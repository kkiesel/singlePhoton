#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
#TODO: look up json ntupelized data and update ANA_JSON
ANA_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt
PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_latest.txt

PU_CALC=pileupCalc.py --calcMode true --minBiasXsec 69300 --maxPileupBin 60 --numPileupBins 60
MAKE_HIST=python makeMChisto.py

puWeights.root: PU_dist.root makeMChisto.py
	$(MAKE_HIST) $@

PU_dist.root: $(ANA_JSON) $(PU_JSON)
	$(PU_CALC) -i $(ANA_JSON) --inputLumiJSON $(PU_JSON) $@


