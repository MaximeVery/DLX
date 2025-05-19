# DLX
Implementation of Knuth's DLX algorithm applied to the search for spreads and packings in finite projective geometry

The file DLX.c is an implementation of the DLX algorithm applied to the computation of spreads and packings in a finite projective space.
It takes the following arguments:

The input file (containing the incidence relations of the projective space)

The output file (which will contain the set of solutions)

The number of threads to execute (for parallelization)

The input file must follow this format: each line corresponds to the labels of the points with which a line is incident (each line of the PG corresponds to a line in the file).

To compile:

gcc -o DLX DLX.c -lm

To run: 

./DLX inputFile.txt outFile.txt numberOfThreads
