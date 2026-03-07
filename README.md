**spacer_count** is a optimized counter for counter spacers or barcodes from sequencing outputs from a variety of methods. It uses regex and pairwise alignment for detecting, therefore it can reliably count spacers appearing at any position, or orientation of the read, even when there are the spacer contains certain level of error (both substitution and in-del). Thus it can be run directly on outputs from sources like whole-plasmid sequencing (long-reads), staggered amplicon, or simple reads without trimming.  

Performance optimization:
1. Multiprocessing if enabled for pairwise alignment. 
2. LRU caching is enabled to accelerate alignment