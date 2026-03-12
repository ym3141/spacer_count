import gzip
import re
import warnings
import time
from multiprocessing import Pool
from functools import lru_cache, partial

import pandas as pd
from Bio import SeqIO, Align, Seq

from os import listdir

class SpacerCounter:

    def __init__(self, flanking_seqs, spacer_size=20, spacer_df=None, spacer_info_csv=None, spacer_size_flex=1):
        left_flanking_seq, right_flanking_seq = flanking_seqs

        if len(left_flanking_seq) < 5 or len(right_flanking_seq) < 5:
            raise Exception("Flanking sequences must be at least 5 bases long.")
        left_flanking_seq = left_flanking_seq.upper()
        right_flanking_seq = right_flanking_seq.upper()
        
        # Warn if there are too many N bases in the flanking sequences, which might lead to too many false positives.
        all_flanking = left_flanking_seq + right_flanking_seq
        if len(all_flanking) - all_flanking.count('N') < 8:
            warnings.warn("Flanking sequences might contain too much N bases ({} Ns).".format(all_flanking.count('N')), UserWarning)
        
        if spacer_info_csv is not None:
            print("Using provided spacer info CSV; ignoring spacer_size parameter.")
            self.spacer_df = pd.read_csv(spacer_info_csv, header=None, names=["guide_id", "sequence", 'gene'])
            self.spacer_df['sequence'] = self.spacer_df['sequence'].str.upper()
            self.spacer_df['guide_id'] = self.spacer_df['guide_id'].astype(str)
            self.spacer_df['gene'] = self.spacer_df['gene'].astype(str)

            self.spacer_size_lims= [self.spacer_df['sequence'].apply(len).min(), self.spacer_df['sequence'].apply(len).max()]
            self.spacer_size_lims = [self.spacer_size_lims[0] - spacer_size_flex, self.spacer_size_lims[1] + spacer_size_flex]

        elif spacer_df is not None:
            print("Using provided spacer info DataFrame; ignoring spacer_size parameter.")
            self.spacer_df = spacer_df
            self.spacer_size_lims= [self.spacer_df['sequence'].apply(len).min(), self.spacer_df['sequence'].apply(len).max()]
            self.spacer_size_lims = [self.spacer_size_lims[0] - spacer_size_flex, self.spacer_size_lims[1] + spacer_size_flex]

        else:
            print("No spacer info DataFrame provided; using spacer_size parameter.")
            self.spacer_size_lims = [spacer_size - spacer_size_flex, spacer_size + spacer_size_flex]
            self.spacer_df = pd.DataFrame(columns=["guide_id", "sequence", 'gene'])


        self.spacer_size_flex = spacer_size_flex

        left_flanking_seq = left_flanking_seq.replace('N', '[ACGT]')
        right_flanking_seq = right_flanking_seq.replace('N', '[ACGT]')

        self.re_pattern = re.compile("{0}((A|C|T|G){{{1},{2}}}){3}".format(
            left_flanking_seq, self.spacer_size_lims[0], self.spacer_size_lims[1], right_flanking_seq))

    '''
    Instance from fasta and csv files. 
    The fasta file should contain two sequences with ids "flanking_left" and "flanking_right". 
    The csv file should contain three columns: guide_id, sequence, and gene (with no header).
    '''
    @classmethod
    def from_fasta_csv(cls, flanking_fasta_path, spacer_info_csv, spacer_size_flex=1):
        with open(flanking_fasta_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id == "flanking_left":
                    left_flanking_seq = str(record.seq).upper()
                elif record.id == "flanking_right":
                    right_flanking_seq = str(record.seq).upper()

        if len(left_flanking_seq) < 5 or len(right_flanking_seq) < 5:
            raise Exception("Flanking sequences must be at least 5 bases long.")

        spacer_df = pd.read_csv(spacer_info_csv, header=None, names=["guide_id", "sequence", 'gene'])
        spacer_df['sequence'] = spacer_df['sequence'].str.upper()
        spacer_df['guide_id'] = spacer_df['guide_id'].astype(str)
        spacer_df['gene'] = spacer_df['gene'].astype(str)

        return cls([left_flanking_seq, right_flanking_seq], 0, spacer_info_csv=spacer_info_csv, spacer_size_flex=spacer_size_flex)


    def count_spacers(self, fastq_path, basename=None, threads=1, first_n=None):
        start_time = time.time()

        # Extract spacers from the fastq file
        id_spacers = self.parse_fastq(fastq_path)
        id_spacers = [(id, spacer) for id, spacer in id_spacers if spacer != ""]

        if first_n is not None:
            id_spacers = id_spacers[:first_n]
            print('  Only processing the first {0} spacers.'.format(first_n))

        print("({1:.2f} seconds)".format(fastq_path, time.time() - start_time))
        start_time = time.time()
        
        # Initialize dict based on spacer_df with an additional row for unknown spacers
        seq_count_dict = {}
        for spacer_seq in self.spacer_df['sequence']:
            if spacer_seq not in seq_count_dict:
                seq_count_dict[spacer_seq] = 0
        unknown_spacer_list = []

        for id, spacer in id_spacers:
            if spacer in seq_count_dict:
                seq_count_dict[spacer] += 1
            elif spacer != "":
                unknown_spacer_list.append((id, spacer))

        print('  Out of total {0} total spacers, {1} ({2:.2%}) matched exactly to a known spacer.'.format(
            len(id_spacers), len(id_spacers) - len(unknown_spacer_list), (len(id_spacers) - len(unknown_spacer_list)) / len(id_spacers)), end=' ')
        print('({1:.2f} seconds)'.format(fastq_path, time.time() - start_time))
        start_time = time.time()

        # sort the reference spacers by their count in descending order, and convert to tuple for caching
        # This way, the most common spacers will be aligned first, reducing the times of alignment needed.
        ref_spaer_list = [k for k, v in sorted(seq_count_dict.items(), key=lambda item: item[1], reverse=True)]
        spacer_tup = tuple(ref_spaer_list)
        align2correct_partial = partial(align2correct, spacer_tup)

        if threads > 1:
            # Use multiprocessing to align unknown spacers in parallel
            with Pool(threads) as pool:
                corrected_spacers = pool.map(align2correct_partial, [spacer for _, spacer in unknown_spacer_list])
            corrected_results = list(zip([id for id, _ in unknown_spacer_list], corrected_spacers))

        else:
            # Align unknown spacers sequentially, this will benefit from lru_cache 
            corrected_results = []
            for id, spacer in unknown_spacer_list:
                corrected_spacer = align2correct_partial(spacer)
                corrected_results.append((id, corrected_spacer))
        
        unknown_spacer_list2 = []
        for idx, (id, spacer) in enumerate(corrected_results):
            if spacer is not None:
                seq_count_dict[spacer] += 1
            else:
                unknown_spacer_list2.append(unknown_spacer_list[idx])

        unknown_dict = {}
        for unknown_id, unknown_seq in unknown_spacer_list2:
            if unknown_seq in unknown_dict:
                unknown_dict[unknown_seq] += 1
            else:
                unknown_dict[unknown_seq] = 1

        
        print('  Out of {1} unknown spacers, {0} were recovered by pairwise alignment.'.format(
            len(unknown_spacer_list) - len(unknown_spacer_list2), len(unknown_spacer_list)), end=' ')
        print('  ({1:.2f} seconds)'.format(fastq_path, time.time() - start_time))

        print('Summary: Out of total {0} spacers, {1} ({2:.2%}) were matched to a known spacer.'.format(
            len(id_spacers), len(id_spacers) - len(unknown_spacer_list2), 
            (len(id_spacers) - len(unknown_spacer_list2)) / len(id_spacers), self.spacer_size_flex))
        print('-' * 40 + '\n')

        output_df = self.spacer_df.copy()
        output_df['count'] = output_df['sequence'].map(seq_count_dict).fillna(0).astype(int)
        output_df.loc[len(output_df.index)] = ['unknown_spacer', 'N' * 5 + '...' + 'N' * 5, 'unknown_gene', len(unknown_spacer_list2)]

        unknown_df = pd.DataFrame(columns=output_df.columns)
        for unknown_seq, count in unknown_dict.items():
            unknown_df = pd.concat([unknown_df, pd.DataFrame([['unknown_spacer', unknown_seq, 'unknown_gene', count]], columns=output_df.columns)], ignore_index=True)

        if basename is not None:
            output_df.to_csv(basename + "spacer_count.csv", index=False)
            unknown_df.to_csv(basename + "unknown_spacer.csv", index=False)

        return output_df, unknown_df

    def parse_fastq(self, fastq_path):

        print("Extracting spacers from file: {0}".format(fastq_path))

        id_spacers = []
        no_guide_count = 0
        id_sequences = load_fasta_to_seqs(fastq_path)
        for id, seq in id_sequences:

            # Try the forward strand
            match = self.re_pattern.search(seq)
            if match:
                id_spacers.append((id, match.group(1)))
                continue

            # Try the reverse strand
            rev_seq = str(Seq.Seq(seq).reverse_complement()).upper()
            match = self.re_pattern.search(rev_seq)
            if match:
                id_spacers.append((id, match.group(1)))
                continue

            # If no match is found, increment the no_guide_count
            id_spacers.append((id, ""))
            no_guide_count += 1

        print('  Out of total {0} reads, {1} ({2:.2%}) likely contain a spacer. (Flexibility = {3})'.format(
            len(id_spacers), len(id_spacers) - no_guide_count, (len(id_spacers) - no_guide_count) / len(id_spacers), self.spacer_size_flex), end=' ')

        return id_spacers
    

@lru_cache(maxsize=1024)
def align2correct(spacer_tup, spacer):
    corrected_spacer = None

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.5

    for index, ref_spacer in enumerate(spacer_tup):
        if aligner.score(spacer, ref_spacer) > 0.9 * len(ref_spacer):
            corrected_spacer = ref_spacer
            break
    
    return (corrected_spacer)
        

def load_fasta_to_seqs(fastq_path):
    id_sequences = []
    if fastq_path.endswith(".gz"):
        with gzip.open(fastq_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                id_sequences.append((record.id, str(record.seq).upper()))
    else:
        with open(fastq_path, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                id_sequences.append((record.id, str(record.seq).upper()))
    return id_sequences

if __name__ == "__main__":
    pass