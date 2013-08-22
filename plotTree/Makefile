VERSION = 13
SUFFIX = _V02.$(VERSION)_tree.root
SAVESUFFIX = new

hadd:
	hadd -f slimGJets$(SUFFIX) slimGJets_200_400$(SUFFIX) slimGJets_400_inf$(SUFFIX)
	hadd -f slimQCD$(SUFFIX) slimQCD_250_500$(SUFFIX) slimQCD_500_1000$(SUFFIX) slimQCD_1000_inf$(SUFFIX)
	hadd -f slimAllQCD$(SUFFIX) slimGJets$(SUFFIX) slimQCD$(SUFFIX)

ewk: ewkClosure ewkFakeRate
	@echo Done

ewkClosure: ewkClosure.py $(addsuffix $(SUFFIX), $(addprefix slim, TTJets WJets))
	@python2 ewkClosure.py --input $(filter-out ewkClosure.py, $^) --savePrefix $(SAVESUFFIX)

ewkFakeRate: ewkFakeRate.py $(addsuffix $(SUFFIX), $(addprefix slim, TTJets WJets DY_to_ee))
	@python2 ewkFakeRate.py --input $(filter-out ewkFakeRate.py, $^) --savePrefix $(SAVESUFFIX)

slim%:% splitCandidates.py
	@python2 splitCandidates.py --input $<  --genMatching

clean:
	rm *~
	rm *.pyc
