# MAPSEQ Analysis

This directory contains our MAPSEQ analysis.
It was converted into a git available repo.

# Instructions

The repo will exclude anything in the `data` directory.
To run the analysis, run the `main.m` file from MATLAB.
If the data is not found, or is not in an expected format, `main.m` will complain.
The analysis will take some time; it's single threaded and resamples data a lot.

# Analysis

The analysis run using this directory is;

* Build a cell-type classifier using the BARseq templates of OB injection dataset.
* Calculates conditional probability of OB injection.
* Calculates conditional probability of PC injection.
* Find circuitry correlation between OB->PC/t.r. and PC->t.r.
