# Single photon analysis
Currently I'm working on the Single Photon SUSY analysis for CMS, so most of the
code will be related to this analysis.

The nTuplizer for SUSY Photon analysis will be used. See the
[twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideSusyNtuplesRA3)
for more infos. Once you have the nTuples, make TTrees out of them
using my TreeWriter.

## TreeWriter
The main part of the analysis. Here trees are generated out ouf the susyNtuples.

## plotTree
The plots scripts allow you to make nice plots out of the trees.
