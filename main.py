from spacer_count.SpacerCounter import SpacerCounter
import pw_align

if __name__ == "__main__":
    counter = SpacerCounter(['NGATG', 'ATGTGGTC'], spacer_size_flex=2, spacer_info_csv='data/spacer_info.csv')
    known_df, unknown_df = counter.count_spacers('data/lr_test.fastq', basename=None, threads=16)

    # print(known_df.head())
    # ref_list = ['AGCTAGCTAGTG', 'TGAAGCCAATAG', 'GGGGTTCCG']

    # print(pw_align.correct_seq(ref_list, 'TGAAGCCTAATAG', 2))
