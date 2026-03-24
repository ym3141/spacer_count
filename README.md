**spacer_count** is a versatile counter for spacers or barcodes from whole-plasimid (or similar type) sequencing outputs. It uses regex and pairwise alignment for sapcer detection, therefore it reliably count spacers appearing in **any** position and orientation in the read, and tolerate error (substitution and in-del). Thus, it can be used directly on vairous outputs from whole-plasmid sequencing (long-reads), staggered amplicon, or simply reads without trimming.  

Performance optimization:
1. Multiprocessing if enabled for pairwise alignment. 
2. LRU caching is enabled to accelerate error correction

# Quick start
To install, simply run:
`pip install spacer_count`

## Input format
The csv file needed for defining known spacer in **spacer_count**, is the same as the [sgRNA file used by MAGeCK](https://sourceforge.net/p/mageck/wiki/input/#h-sgrna-library-file). Briefly, a no header, three-column csv file is used as input.

```
s_10007,TGTTCACAGTATAGTTTGCC,CCNA1
s_10008,TTCTCCCTAATTGCTTGCTG,CCNA1
s_10027,ACATGTTGCTTCCCCTTGCA,CCNC
```

If no input files are provided, all the spacers will be counted in the `unknown` output file.

## Command line interface
**spacer_count**'s primary interface is commend line. After installation, the package is accessible via `spacer-count`.

Example: 

`spacer-count --o data/test --spacer-info-csv data/spacer_info.csv --flanking NGATG-ATGTGGTC -t 4 data/*.fastq`

More argument help is accissbile via `spacer-count --help`

```
usage: spacer_count [-h] --flanking FLANKING [-t THREADS] [--first-n FIRST_N] [-o OUTPUT] [--spacer-flex SPACER_FLEX]
                    [--spacer-info-csv SPACER_INFO_CSV] [--spacer-length SPACER_LENGTH]
                    path [path ...]

Count spacers in fastq files

positional arguments:
  path                  File path(s) to analyze. Can be a single file, or a glob pattern (e.g., 'data/*.fastq'). For 
                        multiple files, separate them with space (e.g.,'data/file1.fastq data/file2.fastq').

options:
  -h, --help            show this help message and exit
  --flanking FLANKING   Flanking (left and right) sequence to consider for spacer counting, e.g., 'NGATG-ATGTGGTC'
  -t THREADS, --threads THREADS
                        Number of threads to use for alignment (default: 1)
  --first-n FIRST_N     Only process the first N sequences in each file (default: all, if not specified)
  -o OUTPUT, --output OUTPUT
                        Base name for output files (default: 'count_table')
  --spacer-flex SPACER_FLEX
                        Allowable flexibility in spacer length in the extracting step (default: 1).
  --spacer-info-csv SPACER_INFO_CSV
                        Path to CSV file containing spacer information (default: None)
  --spacer-length SPACER_LENGTH
                        Expected spacer length. Will be ignored if spacer-info-csv is provided (default: None)
```

## Python interface
Access directely from python is also provide for easier post-processing. The output are `pandas.DataFrame`.

Example: 

```
from spacer_count.SpacerCounter import SpacerCounter

if __name__ == "__main__":
    counter = SpacerCounter(['NGATG', 'ATGTGGTC'], spacer_size_flex=1, spacer_info_csv='data/spacer_info.csv')
    known_df, unknown_df = counter.count_spacers('data/lr_test.fastq', basename=None, threads=8)
    # provide basename to enable file output.

    print(known_df.head())
```
