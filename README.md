# ALGORITHMS FOR BIOINFORMATICS | Blanzieri 

###### Zehra Korkusuz | 2022/2023 

### `Python 3.9` Implementation of Smith-Waterman Algorithm  

Smith-Waterman algorithm is a local alignment algorithm which aims at finding the most similar subsequences between different sequences unlike Needleman-Wush which focuses on finding the optimal global alignment which covers the entire sequence. 

[Background](#background)

[Dependencies](#dependencies)

[How to run the algorithm?](#runalgo)

```console
user@macbook:~$ python main.py TGTTACGG GGTTGACTA --match 3 --mismatch 
-3 --penalty -2

user@macbook:~$ python main.py TGTATTAGCCGG GGTCTGTACTA --match 3 --mismatch 
-3 --penalty -2
```

[Functions used in the algorithm](#functions)

![Smith-Waterman Scoring Gif](https://github.com/zehrakorkusuz/local_sequence_alignment/blob/main/assets/waterman.gif)

### Steps of Smith-Waterman Algorithm {#background}

1. Initilize an empty scoring matrix H with m+1 and n+1 size where m, n are the lengths of the sequences
2. Compute the scores in the matrix with a scoring scheme where default values are
    
    `matching_value = 3`
    `mismatch_value = -3`
   `penalty value = -2`

3. Backtrace recursively from the highest value to smallest value through the path which is greater than 0 to find the optimal local alignment


### APPLICATION {#runalgo}

Go to **smith_waterman** directory

```bash
user@macbook:~$ cd smith_waterman
```
which looks like

```bash
├── smith_waterman/
│   ├── config.py
│   ├── functions.py
│   ├── main.py
│   ├── requirements.txt
```

[config.py](../smith_waterman/config.py) files contains the scoring scheme and can be updated easily

[functions.py](../smith_waterman/functions.py) contains the function required in each step of the algorithm

[main.py](../smith_waterman/main.py)

[requirements.txt](../smith_waterman/requirements.txt) shows the dependencies

### Dependencies {#dependencies}

Download dependencies: 

```bash
user@macbook:~$ pip install numpy tabulate
```

### Running the Application: How to align the sequences? 

Let's consider following sequences for local alignment

`Sequence 1 : TGTTACGG`
`Sequence 1 : GGTTGACTA`

Lets's run the [main.py](../smith_waterman/main.py) file with the scoring arguments and sequences.

- Note that changing the scoring scheme might affect the results and the dafult values are [3,-3, 2] for the algorithm.

```console
user@macbook:~$ python main.py TGTTACGG GGTTGACTA --match 3 --mismatch 
-3 --penalty -2
```
which prints the results. 


```console
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
```

### Functions {#functions}

### `1. smith_waterman(sequence_1, sequence_2)` {#functions}
```python
def smith_waterman(sequence_1, sequence_2):
    """
    The function `smith_waterman` takes two sequences as input and returns the aligned sequences and
    alignment statistics using the Smith-Waterman algorithm.
    
    :param sequence_1: The first sequence of characters or elements that you want to align
    :param sequence_2: The first sequence to be aligned
    :return: three values: base_aligned, match_aligned, and statistics.

    Example:
    """
    H = initialize_matrix(sequence_1, sequence_2)
    # compute the indices of the highest score in the matrix
    i, j = np.unravel_index(np.argmax(H), H.shape)
    base_aligned, match_aligned, statistics = traceback(H, i, j, sequence_1, sequence_2)
    return base_aligned, match_aligned, statistics
```

Smith-Waterman function uses `initialize_matrix` and `traceback `functions

### `2. initialize matrix`

```python
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
```

### `3. traceback`

```python

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
    alignments = []
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
        "best_score" : max_score,
        "matches" : matches,
        "mismatches" : mismatches,
        "gaps" : gaps,
        #"length" : length
    } 

    return alignment_base, alignment_match, statistics
```

### References 

- Image/Gif: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
