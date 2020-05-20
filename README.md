# Chemostat_Rath2017

[![Build Status](https://travis-ci.com/josePereiro/Chemostat_Rath2017.jl.svg?branch=master)](https://travis-ci.com/josePereiro/Chemostat_Rath2017.jl)

## Citing

See `CITATION.bib` for the relevant reference(s).

## (Working) Abstract

Evaluate the framework described in "Fernandez-de-Cossio-Diaz, Jorge, and Roberto Mulet. “Maximum Entropy and Population Heterogeneity in Continuous Cell Cultures.” PLOS Computational Biology 15, no. 2 (February 27, 2019): e1006823. https://doi.org/10.1371/journal.pcbi.1006823" using experimental data from "Rath, Alexander. “Characterisation of Cell Growth, Metabolism and Recombinant Protein Production during Transient and Steady State Conditions for the Human Cell Line AGE1.HN-AAT,” 2017. https://pure.mpg.de/pubman/item/item_2508673_4"

## Package layout

This is a small description of the content in each main folder:

- scripts: Contains the scripts that, if executed in the proper order, compute all the results starting only from the content in `data/raw`.

- src: Contains the source code of the ```Julia``` package itself. In the package I group common tools for the scripts. The scripts will import this package, so you must install it in the environment you'll use to execute them.

- data: Contains all non source code package related files. At cloning, only the `data/raw` folder should contains something.