import pytest
from pw_align import correct_seq_mt
from random import choices

def test_correct_seq():
    ref_list = ['TGAAGCCAATAGTGATTTGGA', 'GAACCGGGCCTAAGCTTTCCG', 'GAGCTAAAATTGCAGCTAGTG']
    assert correct_seq_mt(ref_list, ['GAGCTAAAATTGCAGCTAGTG'], 0, 1) == ['GAGCTAAAATTGCAGCTAGTG']
    assert correct_seq_mt(ref_list, ['GAGCTAAAATTGCAGCTAGTA'], 0, 1) == ['']
    assert correct_seq_mt(ref_list, ['GAGCTAAAATTGCAGCTAGTA'], 1, 1) == ['GAGCTAAAATTGCAGCTAGTG']
    assert correct_seq_mt(ref_list, ['GAGCTAAATTGCAGCTAGTA'], 1, 1) == ['']
    assert correct_seq_mt(ref_list, ['GAGCTAAATTGCAGCTAGTA'], 2, 1) == ['GAGCTAAAATTGCAGCTAGTG']

    assert correct_seq_mt(ref_list, ['TGAAGCCTATAGGATTTGGA'], 2, 1) == ['TGAAGCCAATAGTGATTTGGA']
    assert correct_seq_mt(ref_list, ['TGAAGCCTATAGGATTTGGA'], 1, 1) == ['']

def test_correct_seq_auto():
    ref_list = []
    for _ in range(100):
        ref_list.append(''.join(choices('ACGT', k=22)))

