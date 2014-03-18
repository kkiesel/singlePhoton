#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
#TODO: look up json ntupelized data and update ANA_JSON
ANA_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_latest.txt

PU_CALC=pileupCalc.py --calcMode true --maxPileupBin 60 --numPileupBins 60

puWeights.root: PU_dist.root PU_dist_up.root PU_dist_down.root PU_dist_S10.root
	python2 createWeights.py

PU_dist_S10.root: createMChist.py
	python2 $^ $@

PU_dist.root: $(ANA_JSON) $(PU_JSON) Makefile
	$(PU_CALC) --minBiasXsec 69400 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON) $@

PU_dist_up.root: $(ANA_JSON) $(PU_JSON) Makefile
	$(PU_CALC) --minBiasXsec 72870 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON) $@

PU_dist_down.root: $(ANA_JSON) $(PU_JSON) Makefile
	$(PU_CALC) --minBiasXsec 65930 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON) $@




