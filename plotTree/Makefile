VERSION = 17
SUFFIX = _V02.$(VERSION)_tree.root
SAVESUFFIX = new

AFFIX =

hadd:
	hadd -f slim$(AFFIX)GJets$(SUFFIX) slim$(AFFIX)GJets_200_400$(SUFFIX) slim$(AFFIX)GJets_400_inf$(SUFFIX)
	hadd -f slim$(AFFIX)QCD$(SUFFIX) slim$(AFFIX)QCD_250_500$(SUFFIX) slim$(AFFIX)QCD_500_1000$(SUFFIX) slim$(AFFIX)QCD_1000_inf$(SUFFIX)
	hadd -f slim$(AFFIX)AllQCD$(SUFFIX) slim$(AFFIX)GJets$(SUFFIX) slim$(AFFIX)QCD$(SUFFIX)

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
