import pytest
from spacer_count.SpacerCounter import SpacerCounter, load_fasta_to_seqs, align2correct_mp
import re

class TestExtractSpacers:
    def test_loading_fastq(self):
        id_sequences = load_fasta_to_seqs('data/lr_test.fastq')
        assert len(id_sequences) == 10573
        assert id_sequences[0][0] == '0_SKWKJJ_1'

        id_sequences = load_fasta_to_seqs('data/lr_test.fastq.gz')
        assert len(id_sequences) == 10000
        assert id_sequences[0][0] == '0_NNB8S3_1'

    def test_counter(self):

        with pytest.warns(Warning):
            counter = SpacerCounter(['GANNN', 'ACNNN'], spacer_size_flex=1, spacer_info_csv='data/spacer_info.csv')

        counter = SpacerCounter(['NGATG', 'ATGTGGTC'], spacer_size_flex=1, spacer_info_csv='data/spacer_info.csv')
        assert counter.spacer_df is not None
        assert counter.spacer_df.shape == (11983, 3)
        assert counter.re_pattern == re.compile("[ACGT]GATG((A|C|T|G){35,37})ATGTGGTC")

        output_df, unknown_df = counter.count_spacers('data/lr_test.fastq', basename='data/test_', threads=8, first_n=1000)
        assert output_df.shape[0] == 11984
        assert output_df['count'][9] == 1
        assert unknown_df.shape[0] == 22

    def test_lru_caching(self):
        counter = SpacerCounter(['NGATG', 'ATGTGGTC'], spacer_size_flex=1, spacer_info_csv='data/spacer_info.csv')
        output_df, unknown_df = counter.count_spacers('data/lr_test.fastq', basename='data/test_', threads=1, first_n=1000)
        assert align2correct_mp.cache_info().hits > 0

# test the CLI argument parsing and validation
# uv run -m spacer_count.CLI -o data/test --spacer-info-csv data/spacer_info.csv --flanking NGATG-ATGTGGTC -t 4 data/*.fastq