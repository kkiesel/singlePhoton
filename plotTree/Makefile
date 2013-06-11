VERSION = 12
SUFFIX = _V01.$(VERSION)_tree.root
SAVESUFFIX = new

ewk: ewkClosure ewkFakeRate
	@echo Done

ewkClosure: ewkClosure.py $(addsuffix $(SUFFIX), $(addprefix slim, TTJets WJets))
	@python2 ewkClosure.py --input $(filter-out ewkClosure.py, $^) --savePrefix $(SAVESUFFIX)

ewkFakeRate: ewkFakeRate.py $(addsuffix $(SUFFIX), $(addprefix slim, TTJets WJets DY_to_ee))
	@python2 ewkFakeRate.py --input $(filter-out ewkFakeRate.py, $^) --savePrefix $(SAVESUFFIX)

slim%:% splitCandidates.py
	@python2 splitCandidates.py --genMatching --input $<

