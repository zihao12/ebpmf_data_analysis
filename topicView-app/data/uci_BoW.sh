#!/bin/bash

mkdir uci_BoW
mkdir ../output/uci_BoW
cd uci_BoW
wget http://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/readme.txt
wget http://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/docword.kos.txt.gz
wget http://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/vocab.kos.txt
gunzip docword.kos.txt.gz
