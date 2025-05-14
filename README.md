# Repeat-Aware_Substitution_Rate_Estimator

This tool estimates the substitution rate between two sequences (e.g., two genomes). Our method is robust to the input sequences with high repetitiveness.


# Description of the tool

We consider the following random **substitution process**, parameterized by a rate $0\leq r \leq 1$. Given a string $s$, the character at each position mutates to one of the three other nucleotides with probability $r/3$ per nucleotide independently. The set of all the distinct $k$-mers of the string $s$ is called a **$k$-spectrum** of $s$. The **abundance histogram** of a string $s$ is the sequence $(a_1, \ldots, a_L)$, where $a_i$ is the number of \kmers in $k$-spectrum that occur $i$ times in $s$.

Then, given the $k$-spectrums of $s$ and $t$ and the abundance histogram of $s$, this tool estimates the mutation rate $r$.


# Requirements

Linux (64 bit)

C++17

# installation

```
git clone git@github.com:medvedevgroup/Mutation_rate_estimator.git
cd ./src
make
```

The compiled executable will be located in ./Mutation_rate_estimator/. You can verify successful installation with:

```
./Mutation_rate_estimator -h
```

# Quick Start

We provide example data to help you get started. Run the following commands to test the tool:

## Sequence Mode

```
./Mutation_rate_estimator \
  --mode sequence\
  --input1 ./example_data/origin_seq.fasta \
  --input2 ./example_data/mutated_seq.fasta \
  --k 30
```

## Mixture Mode

Note all the kmers in `./example_data/` are 30-mers.

```
./Mutation_rate_estimator \
  --mode mixture\
  --input1 ./example_data/origin_seq.fasta \
  --input2 ./example_data/mutated_kmers.fasta \
  --k 30
```

## K-mer Mode

```
./Mutation_rate_estimator \
  --mode kmer \
  --input1 ./example_data/origin_kmers.fasta \
  --input2 ./example_data/mutated_kmers.fasta \
  --dist ./example_data/dist.csv \
  --k 30
```

# Usage


## Inputs

### Sequence mode

Provide two complete sequences in separate FASTA files:

```
./Mutation_rate_estimator \
  --mode sequence\
  --input1 seq1.fasta \
  --input2 seq2.fasta \
  --k 31
```


### Mixture mode

User provides a complete sequence and a set of $k$-mers in in separate FASTA files. 

```
./Mutation_rate_estimator \
  --mode mixture\
  --input1 seq1.fasta \
  --input2 set2.fasta \
  --k 31
```

### $K$-mer mode

User provides two sets of $k$-mers in separate FASTA files and an abundance histogram .csv file, in the following form.

```
occurrence i, number of kmers with occurrence i
```

Then use command in folloing form to run

```
./Mutation_rate_estimator \
  --mode kmer \
  --input1 set1.fasta \
  --input2 set2.fasta \
  --dist occ.csv \
  --k 31
```

## Using sketching

User can use parameter `--theta` to speed up the calculation of intersection size of two $k$-spectrums by sketching. We will use estimate the intersection size of the $k$-spectrums of $s$ and $t$ by sketching with sampling rate $\theta$.

The default value of `theta` is $1$, i.e., we select all the $k$-mers to calculate exact intersection size. Note `--theta` should be provided a value in the range of $(0,1]$.

## Other Parameters

`--e`: Absolute error tolerance for Newton’s method, default as $1 \times 10^{-5}$. 

