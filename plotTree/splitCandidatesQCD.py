#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from sys import stdout
from treeFunctions import *

def gammaSelectionClone( tree, newTreeName, type_="tree::Photon", name="photons" ):
	"""Clones a tree and get access to photon vector"""
	newTree = tree.CloneTree(0)
	newTree.SetName( newTreeName )
	photons = ROOT.std.vector(type_)()
	newTree.SetBranchAddress( name, photons )
	# Set the size of the tree in memory to 10MB or get memory overflow.
	newTree.SetMaxVirtualSize(int(1e7))
	return newTree, photons

def getHt( jets, photons, photonJets ):
	ht = 0
	for jet in jets:
		ht += jet.pt
	for photon in photons:
		if not photon._ptJet:
			ht += photon.pt
	for photonJet in photonJets:
		if not photonJet._ptJet:
			ht += photonJet.pt
	return ht

def getMHt( jets, photons, photonJets ):
	mht = ROOT.TVector3(0,0,0)
	j = ROOT.TVector3()
	for jet in jets:
		if not nearPhoton( jet, photonJets ) and not nearPhoton( jet, photons):
			j.SetPtEtaPhi( jet.pt, jet.eta, jet.phi )
			mht += j
	return mht.Pt()

def getMHt2( jets, photons, photonJets ):
	mht = ROOT.TVector3(0,0,0)
	j = ROOT.TVector3()
	for jet in jets:
		if not nearPhoton( jet, photonJets ) and not nearPhoton( jet, photons):
			j.SetPtEtaPhi( jet.pt, jet.eta, jet.phi )
			mht += j
	for jet in photons:
		j.SetPtEtaPhi( jet.ptJet(), jet.eta, jet.phi )
		mht += j
	for jet in photonJets:
		j.SetPtEtaPhi( jet.ptJet(), jet.eta, jet.phi )
		mht += j
	return mht.Pt()

def getHt( jets, photons, photonJets ):
	ht = 0
	for jet in jets:
		ht += jet.pt
	for photon in photons:
		if not photon._ptJet:
			ht += photon.pt
	for photonJet in photonJets:
		if not photonJet._ptJet:
			ht += photonJet.pt
	return ht

def nearPhoton( jet, vector ):
	j = ROOT.TVector3()
	a = ROOT.TVector3()
	j.SetPtEtaPhi( 1, jet.eta, jet.phi )
	for e in vector:
		a.SetPtEtaPhi( 1, e.eta, e.phi )
		if j.DeltaR( a ) < .3:
			return True
	return False

def histoDefinitions():
	hists = {}
	nBins = 3
	hists["eventSplitting"] = ROOT.TH3I("", ";#gamma;#gamma_{jet};#gamma_{e}", nBins, -.5, -.5+nBins, nBins, -.5, -.5+nBins, nBins, -.5, -.5+nBins )
	hists["metSigma"] = ROOT.TH2F("", ";#slash{E}_{T};#sigma_{i#etai#eta}", 50, 0, 500, 440, 0, 0.022 )
	hists["metChIso"] = ROOT.TH2F("", ";#slash{E}_{T};Iso^{#pm}",           50, 0, 500, 300, 0, 30 )
	hists["metNeIso"] = ROOT.TH2F("", ";#slash{E}_{T};Iso^{0}-.04p_{T}",             50, 0, 500, 300, 0, 30 )
	hists["metPhIso"] = ROOT.TH2F("", ";#slash{E}_{T};Iso^{#gamma}-.005p_{T}",        50, 0, 500, 300, 0, 30 )
	hists["metHE"] = ROOT.TH2F("", ";#slash{E}_{T};H/E",                    50, 0, 500, 100, 0, .5 )

	for name, hist in hists.iteritems():
		hist.SetName( name )
		hist.Sumw2()

	return hists


def splitCandidates( inputFileName, processNEvents ):
	"""Key function of splitCanidates.py. The main loop and object selection is
	defined here."""
	datasetAbbr = getDatasetAbbr( inputFileName, slim=False )
	print "Processing file {} with {} configuration".format(inputFileName, datasetAbbr)

	eventHisto = readHisto( inputFileName )
	if processNEvents < 0:
		processNEvents = eventHisto.GetBinContent(1)
	lumiWeight = getLumiWeight( datasetAbbr, processNEvents )

	# Input tree definiton
	tree = readTree( inputFileName, "photonTree" )
	tree.SetName("susyTree") # to avoid overlap with "photonTree" later on

	# Output tree definition
	import os
	outputFileName = "slim"+os.path.basename( inputFileName )
	fout = ROOT.TFile( outputFileName, "recreate" )
	fout.cd()

	photonTree, photons = gammaSelectionClone( tree, "photonTree" )
	photonJetTree, photonJets = gammaSelectionClone( tree, "photonJetTree" )
	photonElectronTree, photonElectrons = gammaSelectionClone( tree, "photonElectronTree" )

	# variables which will be changed in addition to photons:
	ROOT.gROOT.ProcessLine("struct variablesToChange { Float_t ht(0); Float_t weight(0); Float_t st30(0); Float_t st80(0); };")
	change = ROOT.variablesToChange()
	jets = ROOT.std.vector("tree::Jet")()

	for tree_ in [photonTree, photonJetTree, photonElectronTree]:
		tree_.SetBranchAddress("weight", ROOT.AddressOf( change, "weight" ) )
		tree_.SetBranchAddress("ht", ROOT.AddressOf( change, "ht" ) )
		tree_.SetBranchAddress("st30", ROOT.AddressOf( change, "st30" ) )
		tree_.SetBranchAddress("st80", ROOT.AddressOf( change, "st80" ) )
		tree_.SetBranchAddress("jets", jets )

	hists = histoDefinitions()

	for event in tree:
		if not event.GetReadEntry()%100000:
			stdout.write('\r{0} %'.format(100*event.GetReadEntry()/event.GetEntries()))
			stdout.flush()
		if event.GetReadEntry() > processNEvents:
			break

		change.weight = event.weight * lumiWeight
		photons.clear()
		photonElectrons.clear()
		photonJets.clear()
		jets.clear()

		for photon in event.photons:

			# spike rejection
			if photon.r9 <= 1 and photon.sigmaIetaIeta > 0.001:
				# filling the d2 histograms
				if not photon.pixelseed:
					if photon.hadTowOverEm < 0.05 and photon.sigmaIetaIeta < 0.012 and photon.chargedIso < 2.6 and photon.neutralIso < 3.5 +.04*photon.pt:
						hists["metPhIso"].Fill( event.met, photon.photonIso-0.005*photon.pt, change.weight )
					if photon.hadTowOverEm < 0.05 and photon.sigmaIetaIeta < 0.012 and photon.chargedIso < 2.6 and photon.photonIso < 1.3+0.005*photon.pt:
						hists["metNeIso"].Fill( event.met, photon.neutralIso-0.04*photon.pt, change.weight )
					if photon.hadTowOverEm < 0.05 and photon.sigmaIetaIeta < 0.012 and photon.photonIso < 1.3 +0.005*photon.pt and photon.neutralIso < 3.5 +.04*photon.pt:
						hists["metChIso"].Fill( event.met, photon.chargedIso, change.weight )
					if photon.hadTowOverEm < 0.05 and photon.chargedIso < 2.6 and photon.photonIso < 1.3 +0.005*photon.pt and photon.neutralIso < 3.5 +.04*photon.pt:
						hists["metSigma"].Fill( event.met, photon.sigmaIetaIeta, change.weight )
					if photon.sigmaIetaIeta < 0.012 and photon.chargedIso < 2.6 and photon.photonIso < 1.3 +0.005*photon.pt and photon.neutralIso < 3.5 +.04*photon.pt:
						hists["metHE"].Fill( event.met, photon.hadTowOverEm, change.weight )

				# sorting objets
				if photon.hadTowOverEm < 0.05 \
				and photon.sigmaIetaIeta < 0.012 \
				and photon.chargedIso < 2.6 \
				and photon.neutralIso < 3.5 + 0.04*photon.pt \
				and photon.photonIso < 1.3 + 0.005*photon.pt:
					if photon.pixelseed:
						photonElectrons.push_back( photon )
					else:
						photons.push_back( photon )
				elif photon.hadTowOverEm < 0.05 \
				and photon.sigmaIetaIeta < 0.012 \
				and photon.chargedIso < 15 \
				and photon.neutralIso < 20 + 0.04*photon.pt \
				and photon.photonIso < 20 + 0.005*photon.pt \
				and not photon.pixelseed:
					photonJets.push_back( photon )

		change.ht = getHt( event.jets, photons, photonJets )
		if change.ht < 500:
			continue

		for jet in event.jets:
			if not nearPhoton( jet, photons ) and not nearPhoton( jet, photonJets ):
				jets.push_back( jet )
		if jets.size() < 2:
			continue

		change.st30 = getMHt( event.jets, photons, photonJets )
		change.st80 = getMHt2( event.jets, photons, photonJets )

		hists["eventSplitting"].Fill( photons.size(), photonJets.size(), photonElectrons.size(), change.weight )
		if photons.size() and not photonJets.size():
			photonTree.Fill()
		if photonJets.size() and not photons.size():
			photonJetTree.Fill()
		if photonJets.size() and photons.size():
			if photons.at(0).pt > photonJets.at(0).pt:
				photonTree.Fill()
			else:
				photonJetTree.Fill()
		if photonElectrons.size():
			photonElectronTree.Fill()
	print

	# write everything to output file
	photonTree.Write()
	photonJetTree.Write()
	photonElectronTree.Write()

	for name, hist in hists.iteritems():
		hist.Write()

	fout.Close()

if __name__ == "__main__":

	arguments = argparse.ArgumentParser( description="Slim tree" )
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--test", action="store_true" )
	opts = arguments.parse_args()

	# set limit for number of events for testing reason
	processNEvents = 10000 if opts.test else -1

	for inName in opts.filenames:
		splitCandidates( inName, processNEvents )

