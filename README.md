**spacer_count** is a optimized counter for spacers or barcodes from various type of sequencing outputs. It uses regex and pairwise alignment for sapcer detection, therefore it reliably count spacers appearing in **any** position and orientation of the read, and tolerate certain level of error (substitution and in-del). Thus, it can be used directly on outputs from sources like whole-plasmid sequencing (long-reads), staggered amplicon, or simply reads without adaptor trimming.  

Performance optimization:
1. Multiprocessing if enabled for pairwise alignment. 
2. LRU caching is enabled to accelerate error correction

# Quick start
