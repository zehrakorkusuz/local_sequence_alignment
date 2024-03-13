from src.functions import * # See src/functions.py for matrix initialization and traceback
from src.config import * # You can update the scoring scheme from the src/config file
import argparse
import numpy as np
import pandas as pd

"""
    Run the following command in the folder which contains config.py, functions.py and main.py files

    > python main.py TGTTACGG GGTTGACTA --match 3 --mismatch -3 --penalty -2

    SMITH-WATERMAN ALGORITHM FOR LOCAL ALIGNMENT  

    -------------SEQUENCES-------------

    SEQUENCE 1 :  TGTTACGG
    SEQUENCE 2 :  GGTTGACTA 

    -------------ALIGNMENT-------------

    GTT-AC
    ||||||
    GTTGAC
    -------------STATISTICS------------- 

    +------------+---------+------+
    | MISMATCHES | MATCHES | GAPS |
    +------------+---------+------+
    |     0      |    5    |  1   |
    +------------+---------+------+ 

"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Local Alignment with Smith Waterman Algorithm | You can also use config.py file to configure the scoring scheme and filtering values. ')
    

    # Read input sequences from user
    parser.add_argument('sequence_1', type=str, help='Base Sequence')
    parser.add_argument('sequence_2', type=str, help='Sequence to Align with Base Sequence')

    # User can update the scoring scheme
    parser.add_argument('--match', type=int, help='Increment if match is found')
    parser.add_argument('--mismatch', type=int, help='Decrement if a mismatch is found')
    parser.add_argument('--penalty', type=int, help='Penalty value')

    # Filtering options for user
    parser.add_argument('--min_score', type=int, help='Defaul is >=60% of the maximum alignment score')
    parser.add_argument('--min_consecutive_match', type=int, help='At least consecutive match to be considered as local alignment')

    args = parser.parse_args()


    # Store the sequences
    seq1 = args.sequence_1
    seq2 = args.sequence_2

    # Overwrite the scores imported from the config.py if user gives a new scoring scheme
    if args.match is not None:
        scoring['match'] = args.match
    if args.mismatch is not None:
        scoring['mismatch'] = args.mismatch
    if args.penalty is not None:
        scoring['penalty'] = args.penalty
   
    # Overwrite filter config.py if user enters new threshold values on terminal
    if args.min_score is not None:
        filter["min_score"] = args.min_score
    if args.min_consecutive_match is not None:
        filter["min_consecutive_match"] = args.max_min_consecutive_match

    # Run the Smith Waterman Algorithm (See functions.py for more information)

    alignments = smith_waterman(seq1, seq2)


    ### DISPLAY THE ALIGNMENTS

    df = alignment_to_df(alignments)

    ### DISPLAY THE RESULTS ###

    print("\n", "\033[1;97;50m SMITH-WATERMAN ALGORITHM FOR LOCAL ALIGNMENT ", "\n" )

    print( "-------------SEQUENCES-------------")
    print("\n SEQUENCE 1 : ", seq1)
    print(" SEQUENCE 2 : ", seq2, "\n")

    print("-------------ALIGNMENT-------------")
    print(df)