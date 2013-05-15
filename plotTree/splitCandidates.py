#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from treeFunctions import *
import Styles
Styles.tdrStyle2D()
import numpy
import ROOT
import argparse
import ConfigParser
import sys

def deltaPhi( phi1, phi2):
	"""Computes delta phi.
	Due to the zylindrical form of CMS, it has be be corrected by a diffenence
	of 2pi.
	"""
	import math
	result = phi1 - phi2
	while result > math.pi:
		result -= 2*math.pi
	while result <= -math.pi:
		result += 2*math.pi
	return result

def deltaR( object1, object2 ):
	"""Compute delta R of two objects.
	Both objects must have eta and phi as member variables.
	"""
	import math
	return math.sqrt( (object1.eta-object2.eta)**2 + (deltaPhi(object1.phi,object2.phi))**2 )

def gammaSelectionClone( tree, newTreeName, weight, type_="tree::Photon", name="photon" ):
	"""Clones a tree and get access to photon vector"""
	newTree = tree.CloneTree(0)
	newTree.SetName( newTreeName )
	photons = ROOT.std.vector(type_)()
	newTree.SetBranchAddress( name, photons )
	newTree.SetBranchAddress("weight", weight )
	return newTree, photons, weight

def histoDefinition():
	"""All histograms are defined here and put into an directory."""
	from random import randint
	histos = {}
	histos["dEtadPhiLarge"] = ROOT.TH2F( "%x"%randint(0,100000), ";#Delta#eta;#Delta#phi", 100, -4, 4, 100, -3.3, 3.3 )
	histos["dEtadPhiSmall"] = ROOT.TH2F( "%x"%randint(0,1000000), ";|#Delta#eta|;|#Delta#phi|", 100, 0, .012, 100, 0, .1 )
	histos["dEtadPhiSmallPhoton"] = ROOT.TH2F( "%x"%randint(0,1000000), ";|#Delta#eta|;|#Delta#phi|", 100, 0, .012, 100, 0, .1 )
	histos["dEtadPhiSmallElectron"] = ROOT.TH2F( "%x"%randint(0,1000000), ";|#Delta#eta|;|#Delta#phi|", 100, 0, .012, 100, 0, .1 )
	histos["nGenEnRec"] = ROOT.TH2F( "%x"%randint(0,10000000), ";generated e;reco EM-Objects", 5, -0.5, 4.5, 5, -0.5, 4.5 )
	histos["dPtdR"] = ROOT.TH2F( "%x"%randint(0,10000000), ";|(p_{T}-p_{T gen} )/ p_{T gen}|;#DeltaR", 100, 0, 1, 100, 0, .5 )
	return histos

def draw_histogram_dict( histograms, suffix="new" ):
	"""Draw all histograms in this directory and save it as pdf."""
	can = ROOT.TCanvas()
	can.cd()
	dataset = ROOT.TPaveText(.4,.9,.6,.98, "ndc")
	dataset.SetFillColor(0)
	dataset.SetBorderSize(0)
	dataset.AddText( suffix )
	for name, histo in histograms.iteritems():
		name = name.replace(".","")
		histo.Draw("colz")
		dataset.Draw()
		can.SaveAs("plots/%s_%s.pdf"%(name,suffix))

def generalMatching( objects1, objects2, hist):
	"""Matches for each object in objects1 an object in objects2.
	'hists' is a dict containing ROOT.TH* objects which will be filled and returned.
	The gen information will be storend in first list and returned.
	"""
	hist["nGenEnRec"].Fill( objects1.size(), objects2.size() )
	for o1 in objects1:
		match = False
		minDeltaPt = sys.maxint
		minIndex = -1
		for i, o2 in enumerate(objects2):
			DeltaEta = o1.eta - o2.eta
			DeltaPhi = deltaPhi( o1.phi, o2.phi)
			absDeltaEta = abs( DeltaEta )
			absDeltaPhi = abs( DeltaPhi )
			absDeltaPt = 2*abs( o1.pt - o2.pt ) / ( o1.pt + o2.pt )
			hist["dEtadPhiLarge"].Fill( DeltaEta, DeltaPhi )
			hist["dEtadPhiSmall"].Fill( absDeltaEta, absDeltaPhi )
			if hasattr( o1, "pixelseed" ):
				if o1.pixelseed:
					hist["dEtadPhiSmallElectron"].Fill( absDeltaEta, absDeltaPhi )
				else:
					hist["dEtadPhiSmallPhoton"].Fill( absDeltaEta, absDeltaPhi )

			if absDeltaEta < .008 and absDeltaPhi < .08 and absDeltaPt < .1:
				match = True
				if absDeltaPt < minDeltaPt:
					minDeltaPt = absDeltaPt
					minIndex = i

		if match:
			if hasattr( o1, "genInformation" ):
				o1.genInformation = 1
			else: # use phi variable for gen information. TODO: fix that
				# only fill gen matching for photons matching a gen electrons.
				# Electrons are not needed in fake rate
				if not objects2[minIndex].pixelseed:
					o1.phi = 5
	return objects1, hist

def splitCandidates( inputFileName, shortName, nExpected, processNEvents=-1, genMatching=False ):
	print "Processing file {}".format(inputFileName)

	#genMatching = not shortName in ["WJets", "DY", "TTJets"]
	if genMatching:
		print "Match generated objects"

	tree = readTree( inputFileName )
	eventHisto = readHisto( inputFileName )
	if processNEvents == -1:
		processNEvents = eventHisto.GetBinContent(1)

	weight = numpy.zeros(1, dtype=float)
	photonTree, photons, weight = gammaSelectionClone( tree, "photonTree", weight )
	photonJetTree, photonJets, weight = gammaSelectionClone( tree, "photonJetTree", weight )
	photonElectronTree, photonElectrons, weight = gammaSelectionClone( tree, "photonElectronTree", weight )
	if genMatching:
		genElectronTree, genElectrons, weight = gammaSelectionClone( tree, "genElectronTree", weight, "tree::Particle","genElectron" )
		histograms = histoDefinition()

	emObjects = ROOT.std.vector("tree::Photon")()

	for event in tree:
		if not event.GetReadEntry()%100000:
			print '{0}%\r'.format(100*event.GetReadEntry()/event.GetEntries())

		if event.GetReadEntry() > processNEvents:
			break

		weight[0] = event.weight * nExpected / processNEvents
		emObjects.clear()
		photons.clear()
		photonElectrons.clear()
		photonJets.clear()
		if genMatching:
			genElectrons.clear()

		for gamma in event.photon:

			# cuts for every object to rject spikes
			if gamma.r9 < 1 \
			and gamma.sigmaIetaIeta > 0.001 \
			and abs(gamma.eta) < 1.4442 \
			and gamma.pt > 20:

				# look for gamma and electrons
				if gamma.hadTowOverEm < 0.05 \
				and gamma.sigmaIetaIeta < 0.012 \
				and gamma.chargedIso < 2.6 \
				and gamma.neutralIso < 3.5 + 0.04*gamma.pt \
				and gamma.photonIso < 1.3 + 0.005*gamma.pt:
					emObjects.push_back( gamma )

				# QCD fake object definition
				if gamma.ptJet > 75 \
				and gamma.hadTowOverEm < 0.05 \
				and gamma.sigmaIetaIeta < 0.014 \
				and gamma.chargedIso < 15 \
				and gamma.neutralIso < 3.5 + 0.04*gamma.pt \
				and gamma.photonIso < 1.3 + 0.005*gamma.pt \
				and not gamma.pixelseed \
				and ( gamma.sigmaIetaIeta>=0.012 or gamma.chargedIso>=2.6):
					photonJets.push_back( gamma )

		if genMatching:
			emObjects, histograms = generalMatching( emObjects, event.genElectron, histograms )
			for genE in event.genElectron:
				if abs(genE.eta) < 1.4442:
					genElectrons.push_back( genE )
			genElectrons, histograms = generalMatching( genElectrons, emObjects, histograms )

		for emObject in emObjects:
			if emObject.pixelseed:
				photonElectrons.push_back( emObject )
			else:
				photons.push_back( emObject )

		if photons.size() > 0:
			photonTree.Fill()
		else:
			if photonElectrons.size() > 0:
				photonElectronTree.Fill()
			if photonJets.size() > 0:
				photonJetTree.Fill()
		if genMatching and genElectrons.size() > 0:
			genElectronTree.Fill()


	outputFileName = "slim"+inputFileName
	f = ROOT.TFile( outputFileName, "recreate" )
	f.cd()
	photonTree.Write()
	photonJetTree.Write()
	photonElectronTree.Write()
	if genMatching:
		genElectronTree.Write()
		draw_histogram_dict( histograms, shortName )
	f.Close()


if __name__ == "__main__":
	# include knowledge about objects saved in the tree
	ROOT.gSystem.Load("libTreeObjects.so")

	arguments = argparse.ArgumentParser( description="Slim tree" )
	arguments.add_argument("--input", default=["WJets_V01.12_tree.root"], nargs="+" )
	arguments.add_argument("--test", action="store_true" )
	arguments.add_argument("--genMatching", action="store_true" )
	opts = arguments.parse_args()

	processNEvents = 10000 if opts.test else -1


	# disable open canvas
	ROOT.gROOT.SetBatch()

	datasetConfigName = "dataset.cfg"
	datasetConf = ConfigParser.SafeConfigParser()
	datasetConf.read( datasetConfigName )

	integratedLumi = 19300 #pb

	for inName in opts.input:
		shortName = None
		for configName in datasetConf.sections():
			if inName.count( configName ):
				shortName = configName
				crosssection = datasetConf.getfloat( configName, "crosssection" )
		if not shortName:
			print "No configuration for input file {} defined in '{}'".format(
					inName, datasetConfigName )
			continue

		# N = L * sigma
		nExpected = integratedLumi * crosssection

		splitCandidates( inName, shortName, nExpected, processNEvents, opts.genMatching )

