# Entropy Clustering

A C implementation of the Dominance clustering algorithm from [1](#bib1).

## Overview

The following algorithms are implemented for clustering distributions of probabilities with the goal of reducing the partition entropy:

- `rd_clustering`: a simple random assignment of distributions.
- `di_clustering`: The Divisive_Information_Theoretic_Clustering from [2](#bib2).
- `rg_clustering`: The Ratio-Greedy algorithm from [1](#bib1).

The implementations are tested in two data sets: the 20 newsgroup data set available in the scikit-learn package from Python, and the RCV1-v2 data set retrievable from the website http://www.ai.mit.edu/projects/jmlr/papers/volume5/lewis04a/lyrl2004_rcv1v2_README.htm.

## Retrieving data

The data used for the tests presented in Cicalese, Laber, and Murtinho (2019) can be retrieved and saved via the `make_ng20_data.py`and `make_rcv1_data.py` files. The ng20 dataset will be downloaded from scikit-learn, while the files used to prepare the csv files must be downloaded from the RCV1-v2 website:

- http://www.ai.mit.edu/projects/jmlr/papers/volume5/lewis04a/a08-topic-qrels/rcv1-v2.topics.qrels.gz
- http://www.ai.mit.edu/projects/jmlr/papers/volume5/lewis04a/a12-token-files/lyrl2004_tokens_test_pt0.dat.gz
- http://www.ai.mit.edu/projects/jmlr/papers/volume5/lewis04a/a12-token-files/lyrl2004_tokens_test_pt1.dat.gz
- http://www.ai.mit.edu/projects/jmlr/papers/volume5/lewis04a/a12-token-files/lyrl2004_tokens_test_pt2.dat.gz
- http://www.ai.mit.edu/projects/jmlr/papers/volume5/lewis04a/a12-token-files/lyrl2004_tokens_test_pt3.dat.gz
- http://www.ai.mit.edu/projects/jmlr/papers/volume5/lewis04a/a12-token-files/lyrl2004_tokens_train.dat.gz

The files must be downloaded and unzipped in the folder from which `make_rcv1_data.py` is being called.

## Testing

The test file `main_test.c` tests all three algorithms on a single data set for 29 different numbers of clusters (from 2 to 2000). Run from command line as `./test dname n_rows n_cols`:
- `dname`: filename of the csv file with the data set on which the function will be tested (without the `.csv` extension)
- `n_rows`: the humber of rows in the csv file
- `n_cols`: the number of cols in the csv file

If `test` is the name of the executable compiled from new_test.c, then `./test ng20 51480 20` will run the test on the `ng20.csv` file which can be obtained by running `make_ng20_data.py` (see above), while `rcv1 170946 103` will run the test on the provided `rcv1.csv` file. The results will be stored in files `dname_entrs.csv`, `dname_times.csv`, and `dname_iters.csv`, in the folder from which the function is called. The results for the initialization of the DITC algorithm will be stored as well as the results for the first, fifth, tenth, and final iterations of the algorithm.

The test file `single_function_test.c`can be used to manually test a single function in a single data set.

## Bibliography

<a id="#bib1">1</a>: Cicalese, Fernando, Laber, Eduardo, and Murtinho, Lucas. New results on information-theoretic clustering. 2019.

<a id="#bib2">2</a>: Dhillon, Inderjit S., Mallela, Subramanyam, and Kumar, Rahul. A divisive information-theoretic feature clustering algorithm for text classification. *Journal of Machine Learning Research*, 3:1265-1287, 2003.
