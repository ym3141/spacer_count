from spacer_count.SpacerCounter import SpacerCounter


if __name__ == "__main__":
    counter = SpacerCounter(['NGATG', 'ATGTGGTC'], spacer_size_flex=1, spacer_info_csv='data/spacer_info.csv')
    counter.count_spacers('data/lr_test.fastq', basename='data/test_', threads=8)
