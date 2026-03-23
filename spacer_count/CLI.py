import argparse
import sys
from pathlib import Path

from SpacerCounter import SpacerCounter

def main():
    """Main entry point for the SpacerCounter CLI tool."""
    parser = argparse.ArgumentParser(
        description="Count spacers in fastq files",
        prog="spacer_count"
    )
    parser.add_argument("path",
        type=str,
        nargs="+",
        help="File path(s) to analyze. Can be a single file, or a glob pattern (e.g., 'data/*.fastq'). For multiple files, separate them with space (e.g., 'data/file1.fastq data/file2.fastq')."
    )
    parser.add_argument("--flanking",
        required=True,
        type=str,
        help="Flanking (left and right) sequence to consider for spacer counting, e.g., 'NGATG-ATGTGGTC'"
    )
    parser.add_argument("-t", "--threads",
        type=int,
        default=1,
        help="Number of threads to use for alignment (default: 1)"
    )
    parser.add_argument("--first-n",
        type=int,
        default=None,
        help="Only process the first N sequences in each file (default: all, if not specified)"
    )
    parser.add_argument("-o", "--output",
        type=str,
        default="count_table",
        help="Base name for output files (default: 'count_table')"
        )
    parser.add_argument("--spacer-flex",
        type=int,
        default=0,
        help="Allowable flexibility in spacer length in the extracting step (default: 0)."
    )
    parser.add_argument("--spacer-info-csv",
        type=str,
        default=None,
        help="Path to CSV file containing spacer information (default: None)"
    )
    parser.add_argument("--spacer-length",
        type=int,
        default=None,
        help="Expected spacer length. Will be ignored if spacer-info-csv is provided (default: None)"
    )

    
    args = parser.parse_args()

    # validate the input files
    fastq_files = []    
    for path in args.path:
        if not Path(path).is_file():
            print(f"Error: File '{path}' does not exist.", file=sys.stderr)
            sys.exit(1)
        elif path.endswith('.fastq') or path.endswith('.fastq.gz'):
            fastq_files.append(path)

    # validate the flanking argument
    if '-' not in args.flanking:
        print("Error: Flanking sequences must be provided in the format 'LEFT-RIGHT', e.g., 'NGATG-ATGTGGTC'.", file=sys.stderr)
        sys.exit(1)
    else:
        left_flank, right_flank, = args.flanking.split('-')
        if not left_flank or not right_flank:
            print("Error: Both left and right flanking sequences must be provided in the format 'LEFT-RIGHT', e.g., 'NGATG-ATGTGGTC'.", file=sys.stderr)
            sys.exit(1)

    # make sure either spacer_info_csv or spacer_length is provided
    if not args.spacer_info_csv and not args.spacer_length:
        print("Error: Either --spacer-info-csv or --spacer-length must be provided.", file=sys.stderr)
        sys.exit(1)
    
    # Construct the SpacerCounter and run the counting
    counter = SpacerCounter((left_flank, right_flank), spacer_size_flex=args.spacer_flex, spacer_info_csv=args.spacer_info_csv, spacer_len=args.spacer_length)

    for fastq_file in fastq_files:
        counter.count_spacers(fastq_file, basename=args.output, threads=args.threads, first_n=args.first_n)

if __name__ == "__main__":
    main()
