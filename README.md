# Mutation_rate_estimator

The tool to estimate substitution rate based on k-spectrum



# Usage

## Inputs

### Sequence mode

User provides two whole sequences in separated fasta files.

### Mixture mode

User provides a whole sequence and a set of $k$-mers in a fasta file. 

### $K$-mer mode

User provides two sets of  of $k$-mers in a fasta file and an additional distribution csv file, in the following form.

```
occurrence i, number of kmers with occurrence i
```

## Using sketching

Set a threshold $0\leq \theta \leq1$. 
