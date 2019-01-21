#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer
import numpy as np

# fetch_20newsgroups parameters
SUBSET = 'all'
RANDOM_STATE=1
REMOVE = ['headers', 'quotes', 'footers']

# CountVectorizer parameters
MAX_DF = 0.95
MIN_DF = 2
STOP_WORDS = 'english'

# fetch and vectorize data

def make_ng20_file():
    ng20 = fetch_20newsgroups(subset=SUBSET, random_state=RANDOM_STATE,
                              remove=REMOVE)

    cv = CountVectorizer(max_df=MAX_DF, min_df=MIN_DF, stop_words=STOP_WORDS)

    vec_data = cv.fit_transform(ng20.data)
    n_docs, n_words = vec_data.shape
    n_classes = len(np.unique(ng20.target))

    # probability of classes
    docs_per_class = np.unique(ng20.target, return_counts=True)[1]
    prob_class = docs_per_class / docs_per_class.sum()

    # words per class (smoothed)
    nums = np.ones((n_words, n_classes))
    for i in range(n_docs):
        d = vec_data[i]
        c = ng20.target[i]
        for w in d.indices:
            nums[w, c] += d[0,w]

    # probability of classes given word
    prob_words_given_class = nums / nums.sum(axis=0)
    ans = prob_words_given_class * prob_class
    np.savetxt("ng20.csv", ans, delimiter=",")

if __name__ == "__main__":
    make_ng20_file()
