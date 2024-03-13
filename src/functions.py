import numpy as np
import pandas as pd
from src.config import *

def smith_waterman(sequence_1, sequence_2):
    """
    The function `smith_waterman` takes two sequences as input and returns the aligned sequences and
    alignment statistics using the Smith-Waterman algorithm.
    
    :param sequence_1: The first sequence of characters or elements that you want to align
    :param sequence_2: The first sequence to be aligned
    :return: three values: base_aligned, match_aligned, and statistics.

    INPUT:

    seq1 = "TGTTACGG" 
    seq2 = "GGTTGACTA"  

    OUTPUT:

    ('GTT-AC',
    'GTTGAC',
    {'best_score': 13.0, 'matches': 5, 'mismatches': 0, 'gaps': 1})

    """
    H = initialize_matrix(sequence_1, sequence_2)
    alignments = []

    # find indices to traceback 
    indices = find_indices(H, threshold=filter["min_score"]) # which has threshold value extracted from config, if you want you can set by threshold=10 etc. 

    # traceback and add to alignment if meet the requirement given under the config by user
    for i in indices:
        i, j = i[0], i[1]
        result = traceback(H, i, j, sequence_1, sequence_2)
        condition = filter_alignments(result["stats"])
        if condition:
            alignments.append(result)

    return alignments

def initialize_matrix(base_sequence="AGCTA", matching_sequence="GCTAA", match_value=scoring["match"], mismatch_value=scoring["mismatch"], penalty_value=scoring["penalty"]):
    """
    Example:

    Input: 
    seq1 = "TGTTACGG" 
    seq2 = "GGTTGACTA"  

    Output:
    array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  3.,  3.,  1.,  0.,  0.,  3.,  1.],
       [ 0.,  3.,  3.,  1.,  1.,  6.,  4.,  2.,  1.,  0.],
       [ 0.,  1.,  1.,  6.,  4.,  4.,  3.,  1.,  5.,  3.],
       [ 0.,  0.,  0.,  4.,  9.,  7.,  5.,  3.,  4.,  2.],
       [ 0.,  0.,  0.,  2.,  7.,  6., 10.,  8.,  6.,  7.],
       [ 0.,  0.,  0.,  0.,  5.,  4.,  8., 13., 11.,  9.],
       [ 0.,  3.,  3.,  1.,  3.,  8.,  6., 11., 10.,  8.],
       [ 0.,  3.,  6.,  4.,  2.,  6.,  5.,  9.,  8.,  7.]])

    """
    m, n = len(base_sequence), len(matching_sequence)
    H = np.zeros((m+1, n+1))
    for i in range(1, m+1):
        for j in range(1, n+1):
            base_score = H[i-1,j-1]
            #print(base_score)
            left = H[i, j-1] + penalty_value
            up = H[i-1, j] + penalty_value
            if base_sequence[i-1] == matching_sequence[j-1]:
                #print(seq1[i-1], seq2[j-1])
                score = base_score + match_value
            else:
                score = base_score + mismatch_value

            H[i,j] = max(0, score, left, up)
    return H

# traceback, recursive function to traceback to find the local alignment 
def traceback(H, i, j, seq1, seq2):
    """
    Given the location of the highest score, tracebacks till the beginning of the local alignment

    INPUT:

    H = array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  3.,  3.,  1.,  0.,  0.,  3.,  1.],
       [ 0.,  3.,  3.,  1.,  1.,  6.,  4.,  2.,  1.,  0.],
       [ 0.,  1.,  1.,  6.,  4.,  4.,  3.,  1.,  5.,  3.],
       [ 0.,  0.,  0.,  4.,  9.,  7.,  5.,  3.,  4.,  2.],
       [ 0.,  0.,  0.,  2.,  7.,  6., 10.,  8.,  6.,  7.],
       [ 0.,  0.,  0.,  0.,  5.,  4.,  8., 13., 11.,  9.],
       [ 0.,  3.,  3.,  1.,  3.,  8.,  6., 11., 10.,  8.],
       [ 0.,  3.,  6.,  4.,  2.,  6.,  5.,  9.,  8.,  7.]])

    i = 6
    j = 7
    seq1 = "TGTTACGG" 
    seq2 = "GGTTGACTA"  

    OUTPUT: 

    ('GTT-AC',
    'GTTGAC',
    {'best_score': 13.0, 'matches': 5, 'mismatches': 0, 'gaps': 1})

    """
    alignment_base = ""
    alignment_match = ""
    max_score = H[i][j]
    gaps, mismatches, matches = 0, 0, 0
    while H[i][j] != 0:
        if seq1[i-1] == seq2[j-1]:
            #print("match", seq1[i-1], "-" ,seq2[j-1])
            alignment_base = seq1[i-1] + alignment_base
            alignment_match = seq1[i-1] + alignment_match
            i -= 1
            j -= 1
            matches += 1
        else:
            if H[i][j] == H[i-1][j-1] + scoring["mismatch"]:
                alignment_base = seq1[i-1] + alignment_base
                alignment_match = seq2[j-1] + alignment_match
                i -= 1
                j -= 1
                mismatches += 1
            elif H[i][j] == H[i][j-1] + scoring["penalty"]:
                #print("go up up up", )
                alignment_base = "-" + alignment_base
                alignment_match = seq2[j-1] + alignment_match
                j -= 1
                gaps += 1
            else: #mx[i][j] == mx[i-1][j] -2:
                #print("go left")
                alignment_base = seq1[i-1] + alignment_base
                alignment_match = "-" + alignment_match
                i -= 1
                gaps += 1

    statistics = {
        "score" : max_score,
        "matches" : matches,
        "mismatches" : mismatches,
        "gaps" : gaps,
    } 

    return {"seq1": alignment_base, "seq2":alignment_match, "stats": statistics} # indices to be iterable

def find_indices(matrix, threshold=filter["min_score"]):

    high_scoring_indices = np.argwhere(matrix >= threshold)
    high_scoring_values = matrix[high_scoring_indices[:, 0], high_scoring_indices[:, 1]]

    #sort by value of the indices in the matrix
    sorted_indices = high_scoring_indices[np.argsort(-high_scoring_values)]
    sorted_values = high_scoring_values[np.argsort(-high_scoring_values)]

    #filter values by threshold
    above_threshold_mask = sorted_values >= threshold
    sorted_indices = sorted_indices[above_threshold_mask]

    return sorted_indices

def filter_alignments(stats, min_score= filter["min_consecutive_match"], min_consecutive_match=filter["min_consecutive_match"]):
    if min_score is not None:
        if stats["score"] < min_score:
            return False
    if min_consecutive_match is not None:
        if stats["max_consecutive_match"] < min_consecutive_match :
            return False
    return True

def alignment_to_df(alignments):
        df = pd.DataFrame(alignments)
        df_stats = pd.json_normalize(df["stats"])
        df_stats = df_stats.rename(columns={
            "best_score" : "Score",
            "matches" : "Matches",
            "mismatches" : "Mismatches",
            "gaps" : "Gaps" 
        })

        df = pd.concat([df[['seq1', 'seq2']], df_stats], axis=1)
         # Set index starting from 0, to later access easier
        df.index = range(len(df))

        return df