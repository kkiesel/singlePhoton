# TreeWriter
The main part of the analysis.

## Preparations
Here susyTrees are generated out ouf the susyNtuples.
The nTuples are produced by the [nTuplizer](https://github.com/CMSSUSYPhotons/SUSYPhotonAnalysis)
for SUSY Photon analysis. See the [twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideSusyNtuplesRA3)
for more information.

## Using the treeWriter
Compile with make. For running on MC, make sure that the input for the
pileUpReweigting is pointing to some files. Default is on */afs/cern.ch*, so this
should not be a problem.

If you want to do use the susyTrees afterwards, it is useful to copy the
*libTreeObjects.so* to *$PATH* or include this folder in *$PATH*.

The important settings can be made in *executable.cc* (recompile):
*SplitTrees()* will allow you to either store all photons in the photonCollection,
or split the events in tight, loose and electron-like-photon events.

Use
    ./execute outputFileName.root input1.root input2.root ...
to create susyTrees.

## Running on NAF
Comment out the datasets you want to process in *submitAll*, set the version and
execute the script.
