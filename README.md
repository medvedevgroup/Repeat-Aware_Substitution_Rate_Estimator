# Mutation_rate_estimator

This tool estimate substitution rate based on k-spectrum.Let $s$ be a string and let $t$ be generated from $s$ using a substitution process with mutation rate $r$. Given the $k$-spectrums of $s$ and $t$ and the abundance histogram of $s$, this tool estimates the mutation rate $r$.

# Requirements

Linux (64 bit)
C++17

# installaion

```
git clone git@github.com:medvedevgroup/Mutation_rate_estimator.git
cd ./src
make
```

Note the executable file is located in `./Mutation_rate_estimator/`. You can use following command to check if you install the tool successfully. 

```
./Mutation_rate_estimator -h
```

# Usage


## Inputs

### Sequence mode

User provides two whole sequences in separated fasta files.

```
./Mutation_rate_estimator \
  --mode sequence\
  --input1 seq1.fasta \
  --input2 seq2.fasta \
  --k 31
```


### Mixture mode

User provides a whole sequence and a set of $k$-mers in a fasta file. 

```
./Mutation_rate_estimator \
  --mode mixture\
  --input1 seq1.fasta \
  --input2 set2.fasta \
  --k 31
```

### $K$-mer mode

User provides two sets of $k$-mers in a fasta file and an additional distribution csv file, in the following form.

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

User can use parameter `--theta` to speed up the calculation of intersection size of two $k$-spectrums by sketching. The default value of `theta` is $1$, i.e., we select all the $k$-mers to calculate exact intersection size. Note `--theta` should be provided a value in the range of $(0,1]$.

## Other Parameters

1. `--e`: the abosulte error upper bound for Newton's Method, default as $1 \times 10^{-5}$. 

2. TBD ...


# Example

We prepare exmaples to run the code. You can run following commands to test the tool.

```
./Mutation_rate_estimator \
  --mode sequence\
  --input1 ./example_data/origin_seq.fasta \
  --input2 ./example_data/mutated_seq.fasta \
  --k 30
```

Note all the kmers in `./example_data/` are 30-mers.

```
./Mutation_rate_estimator \
  --mode mixture\
  --input1 ./example_data/origin_seq.fasta \
  --input2 ./example_data/mutated_kmers.fasta \
  --k 30
```

```
./Mutation_rate_estimator \
  --mode kmer \
  --input1 ./example_data/origin_kmers.fasta \
  --input2 ./example_data/mutated_kmers.fasta \
  --dist ./example_data/dist.csv \
  --k 30
```
