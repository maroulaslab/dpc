This file provides an implementation in R of the dpc distance for persistence diagrams 
introduced by Marchese and Maroulas in the 2018 paper Signal classification with a point 
process distance on the space of persistence diagrams.

Input:  d1, d2: persistence diagrams, as created by the R package 'TDA'
        beta: a natural number representing the homological dimension
        c: a positive number representing the penalty for cardinality mismatch
        p: a positive number for the exponent used in calculating the distance
        
