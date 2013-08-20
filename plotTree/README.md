# Plotting scripts

## Preparations
Produce susyTrees (*../treeWriter/*) and have a close look at (*../treeWriter/treeWriter.cc*)
Make sure the *libTreeObjects.so* is in *$PATH*.

## splitCandidatesQCDtest.py splitCandidates.py
If the splitting is not done at the treeWriter step, it has to be done here.
HT has to be recalculated and the jets have to be cleaned from the photon(-like)
objects. Additionaly, the crosssections are read from *datasets.cfg* and the
samples are scaled by the luminosity.

You can now add the samples with *hadd*.

## qcdClosureV02.py
QCD closure test

## drawHistosIsoDependency.py
Reads 2d Histos and illustrate the content in 1d plots
