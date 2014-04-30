# Null models of bipartite network structure in C

Timothée Poisot - timothee.poisot@uqar.ca

Requires the C GSL library

Compile with `make`, or just run

~~~ bash
gcc nulls.c -o nm -lgsl -lgslcblas -O3 -DHAVE_INLINE
~~~

Run with ```./nm test.txt 3 7 1000 2000```

## Arguments

1. web file
1. number of columns
1. number of rows
1. target number of replicates
1. maximal number of iterations 
