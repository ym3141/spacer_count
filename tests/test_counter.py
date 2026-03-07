import pytest
from src.spacer_counter import SpacerCounter, load_fasta_to_seqs, align2correct
import re

class TestExtractSpacers:
    def test_loading_fastq(self):
        id_sequences = load_fasta_to_seqs('data/long_read.fastq')
        assert len(id_sequences) == 10000
        assert id_sequences[0][0] == '0_NNB8S3_1'

        id_sequences = load_fasta_to_seqs('data/long_read.fastq.gz')
        assert len(id_sequences) == 10000
        assert id_sequences[0][0] == '0_NNB8S3_1'

    def test_counter(self):
        counter = SpacerCounter(['GATCT', 'ACGCG'], spacer_size_flex=1, spacer_info_csv='data/spacer_info.csv')
        assert counter.spacer_df is not None
        assert counter.spacer_df.shape == (254, 3)
        assert counter.re_pattern == re.compile("GATCT((A|C|T|G){19,21})ACGCG")

        output_df, unknown_df = counter.count_spacers('data/long_read.fastq', basename='data/test_', threads=8)
        assert output_df.shape[0] == 255
        assert output_df['count'][2] == 247
        assert unknown_df.shape[0] == 149

    def test_lru_caching(self):
        counter = SpacerCounter(['GATCT', 'ACGCG'], spacer_size_flex=1, spacer_info_csv='data/spacer_info.csv')
        output_df, unknown_df = counter.count_spacers('data/long_read.fastq', basename='data/test_', threads=1)
        assert align2correct.cache_info().hits > 0
