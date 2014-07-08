#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
ANA_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_latest.txt

PU_CALC=pileupCalc.py --calcMode true --maxPileupBin 60 --numPileupBins 60 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON)

puWeights.root: createWeights.py nTrueVertexMC.root nTrueVertexData.root
	python createWeights.py --mc nTrueVertexMC.root --data nTrueVertexData.root -o $@

nTrueVertexMC.root: createMChist.py
	python $^ $@

# for the pileup uncertainties shift cross section by 5% as proposed in
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVPileUpDescription#Pileup_Systematic_Errors
# and https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupSystematicErrors
nTrueVertexData.root: $(ANA_JSON) $(PU_JSON) Makefile
	$(PU_CALC) --minBiasXsec 69400 --pileupHistName pileup tmpHist1.root
	$(PU_CALC) --minBiasXsec 72870 --pileupHistName pileupUp tmpHist2.root
	$(PU_CALC) --minBiasXsec 65930 --pileupHistName pileupDown tmpHist3.root
	hadd -f $@ tmpHist1.root tmpHist2.root tmpHist3.root
	rm tmpHist1.root tmpHist2.root tmpHist3.root




