#!/usr/bin/env python3
"""
Anderson_deduper.py

Reference-based PCR duplicate remover for single-end reads with known UMIs.

Assumptions:
- Input SAM is coordinate-sorted (use `samtools sort` before running).
- Reads are uniquely mapped (unmapped/secondary/supplementary skipped).
- UMIs are known and provided in a text file (one per line).
- UMI = last ':'-separated token in QNAME.
- Duplicate key = (RNAME, strand, 5'-position, UMI)
  * + strand: 5' = POS - leading_soft_clip
  * - strand: 5' = (POS + ref_consumed(CIGAR) - 1) + trailing_soft_clip
  * ref_consumed counts M, =, X, D, N only.
- First read encountered per key is written.
- Memory-efficient: clears duplicates when chromosome changes.
"""

import argparse
import re

# Constants for FLAG bits
FLAG_REVERSE = 0x10
FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800

# Regex for parsing CIGAR strings
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=XB])")


def get_args():
    """Parse command-line arguments, including a manual -h/--help option."""
    parser = argparse.ArgumentParser(
        prog="Anderson_deduper.py",
        add_help=False
    )

    parser.add_argument("-f", "--file", required=False,
                        help="Path to sorted input SAM file.")
    parser.add_argument("-o", "--outfile", required=False,
                        help="Path to deduplicated output SAM file.")
    parser.add_argument("-u", "--umi", required=False,
                        help="File containing known UMIs (one per line).")
    parser.add_argument("-h", "--help", action="store_true",
                        help="Show this help message and exit.")

    args = parser.parse_args()

    if args.help or not (args.file and args.outfile and args.umi):
        print("""
Usage:
    ./Anderson_deduper.py -u <UMI_file> -f <sorted_input.sam> -o <output.sam>

Description:
    Removes PCR duplicates from a coordinate-sorted SAM file of uniquely mapped
    single-end reads with known UMIs.

Assumptions:
    - The SAM file is coordinate-sorted (use `samtools sort`).
    - Reads are uniquely mapped; unmapped, secondary, and supplementary reads
      are skipped.
    - UMI is the last ':'-separated token in QNAME.
    - Duplicate key: (RNAME, strand, 5'-position, UMI).
    - First read per key is retained.
    - Memory-efficient: clears duplicates when chromosome changes.

Example:
    ./Anderson_deduper.py -u STL96.txt -f input.sorted.sam -o output.dedup.sam

Required Arguments:
    -u, --umi       Path to file containing known UMIs (one per line)
    -f, --file      Path to coordinate-sorted input SAM file
    -o, --outfile   Path to deduplicated output SAM file
    -h, --help      Show this help message and exit
        """)
        exit()

    return args


def load_umis(path):
    """Load whitelist of UMIs from text file."""
    umis = set()
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                umis.add(line)
    return umis


def is_unmapped(flag):
    return flag & FLAG_UNMAPPED != 0


def is_secondary_or_supp(flag):
    return (flag & FLAG_SECONDARY) or (flag & FLAG_SUPPLEMENTARY)


def strand_from_flag(flag):
    return "-" if (flag & FLAG_REVERSE) else "+"


def ref_consumed_from_cigar(cigar):
    """Return number of reference bases consumed by CIGAR string."""
    total = 0
    for length, op in CIGAR_RE.findall(cigar):
        n = int(length)
        if op in ("M", "=", "X", "D", "N"):
            total += n
    return total


def compute_five_prime(flag, pos, cigar):
    """
    Compute 5' genomic coordinate with soft-clip adjustment at the 5' end
    (in read orientation).

    + strand:
        five_prime = POS - leading_soft_clip
    - strand:
        end = POS + ref_consumed(CIGAR) - 1
        five_prime = end + trailing_soft_clip

    Notes:
      - leading_soft_clip = S at the start of CIGAR (e.g., "5S95M...")
      - trailing_soft_clip = S at the end of CIGAR (e.g., "...95M5S")
      - S does not consume reference; we adjust explicitly for true fragment start.
    """
    parts = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    leading_s = parts[0][0] if parts and parts[0][1] == "S" else 0
    trailing_s = parts[-1][0] if parts and parts[-1][1] == "S" else 0

    if flag & FLAG_REVERSE:
        end = pos + ref_consumed_from_cigar(cigar) - 1
        return end + trailing_s
    else:
        return pos - leading_s


def extract_umi_from_qname(qname):
    """Return last ':'-separated token as UMI."""
    parts = qname.strip().split(":")
    return parts[-1]


def write_headers_and_get_first(fin, fout):
    """
    Write SAM headers to output; return (first non-header line or None, header_count).
    """
    first_line = None
    header_count = 0
    for line in fin:
        if line.startswith("@"):
            fout.write(line)
            header_count += 1
        else:
            first_line = line
            break
    return first_line, header_count


def dedupe_stream(in_path, out_path, umi_set):
    """Stream through SAM file and remove PCR duplicates; print formatted summary."""
    kept = dupes = bad_umi = 0
    current_chr = None
    seen = set()
    chr_counts = {}
    header_count = 0

    with open(in_path, "r") as fin, open(out_path, "w") as fout:
        first_line, header_count = write_headers_and_get_first(fin, fout)
        if first_line:
            lines = [first_line]
        else:
            lines = []

        for line in lines + list(fin):
            if line.startswith("@"):
                # Normally headers won't appear here, but be tolerant.
                fout.write(line)
                header_count += 1
                continue

            cols = line.strip().split("\t")
            if len(cols) < 11:
                continue

            qname = cols[0]
            flag = int(cols[1])
            rname = cols[2]
            pos = int(cols[3])
            cigar = cols[5]

            # Skip unmapped / secondary / supplementary
            if is_unmapped(flag) or is_secondary_or_supp(flag):
                continue

            # Reset seen when chromosome changes (sorted input assumed)
            if current_chr is None:
                current_chr = rname
            elif rname != current_chr:
                seen.clear()
                current_chr = rname

            umi = extract_umi_from_qname(qname)
            if umi not in umi_set:
                bad_umi += 1
                continue

            strand = strand_from_flag(flag)
            five_prime = compute_five_prime(flag, pos, cigar)
            key = (rname, strand, five_prime, umi)

            if key in seen:
                dupes += 1
                continue

            # Keep first observation of this key
            fout.write(line if line.endswith("\n") else line + "\n")
            seen.add(key)
            kept += 1
            chr_counts[rname] = chr_counts.get(rname, 0) + 1

    # ---- Summary output in requested format ----
    print("headers\t" + str(header_count))
    print("unique_reads\t" + str(kept))
    print("invalid_umi_records\t" + str(bad_umi))
    print("duplicates_removed\t" + str(dupes))

    # Simple lexicographic order for chromosomes
    for chrom in sorted(chr_counts.keys()):
        print(f"{chrom}\t{chr_counts[chrom]}")


# ---- Run script ----
args = get_args()
umi_set = load_umis(args.umi)
if not umi_set:
    print("Error: UMI list empty or unreadable.")
else:
    dedupe_stream(args.file, args.outfile, umi_set)

