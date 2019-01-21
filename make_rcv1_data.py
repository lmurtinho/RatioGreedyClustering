#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 21:31:44 2019

@author: lucasmurtinho
"""

from sklearn.feature_extraction.text import CountVectorizer
import numpy as np
import random
import time


random.seed(1)

def get_rcv1_vectors(doc_files, class_file):
    t1 = time.time()
    file_dict = {}
    for file in doc_files:
        print("file:", file)
        with open(file, "r") as f:
            s = ""
            for l in f.readlines():
                if l[0] == ".":
                    token_list = l.split(" ")
                    if len(token_list) > 1:
                        key = token_list[1][:-1]
                elif l[0] != "\n":
                    s += l[:-1] + " "
                elif s:
                    file_dict[key] = s
                    s = ""
            if s:
                file_dict[key] = s
        print(len(file_dict), "total elements")
    t2 = time.time()
    print("files to dict in {} seconds".format(t2-t1))


    t1 = time.time()
    class_dict = {key: "" for key in file_dict.keys()}
    with open(class_file, 'r') as f:
        for l in f.readlines():
            val, key, _ = l.split(" ")
            class_dict[key] += val + " "
    t2 = time.time()

    print("classes to dict in {} seconds".format(t2-t1))

    t1 = time.time()
    MAX_DF = 0.95
    MIN_DF = 2
    STOP_WORDS = 'english'
    cv = CountVectorizer(max_df=MAX_DF, min_df = MIN_DF,
                         stop_words=STOP_WORDS)
    data = cv.fit_transform(file_dict.values())
    t2 = time.time()
    print("data vectorization done in {} seconds".format(t2-t1))

    t1 = time.time()
    class_data = cv.fit_transform(class_dict.values())
    t2 = time.time()
    print("class vectorization done in {} seconds".format(t2-t1))

    return data, class_data

def choose_classes(class_data):
    return [random.choice(class_data[i].indices)
            for i in range(class_data.shape[0])]

def prob_class_given_word(data, classes):

    n_docs, n_words = data.shape
    print("{} docs and {} words".format( n_docs, n_words))

    # probability of classes
    docs_per_class = np.unique(classes, return_counts=True)[1]
    prob_class = docs_per_class / docs_per_class.sum()

    n_classes = len(docs_per_class)
    print("{} classes".format(n_classes))

    # words per class (smoothed)
    nums = np.ones((n_words, n_classes))
    for i in range(n_docs):
        if not (i % 100):
            print("iteration", i)
        d = data[i]
        c = classes[i]
        for w in d.indices:
            nums[w, c] += d[0,w]
            #words_per_class[c] += d[0, w]

    # probability of classes given word
    prob_words_given_class = nums / nums.sum(axis=0)
    ans = prob_words_given_class * prob_class
    return ans

if __name__ == "__main__":

    doc_files = ['lyrl2004_tokens_train.dat',
                 'lyrl2004_tokens_test_pt0.dat',
                 'lyrl2004_tokens_test_pt1.dat',
                 'lyrl2004_tokens_test_pt2.dat',
                 'lyrl2004_tokens_test_pt3.dat']

    class_file = 'rcv1-v2.topics.qrels'

    t1 = time.time()
    data, class_data = get_rcv1_vectors(doc_files, class_file)
    t2 = time.time()
    print("Vectorization done in {} seconds".format(t2-t1))

    t1 = time.time()
    classes = choose_classes(class_data)
    t2 = time.time()
    print("Classes randomly chosen in {} seconds".format(t2-t1))

    t1 = time.time()
    ans = prob_class_given_word(data, classes)
    t2 = time.time()
    print("Final array built in {} seconds".format(t2-t1))
    np.savetxt("rcv1.csv", ans, delimiter=",")
