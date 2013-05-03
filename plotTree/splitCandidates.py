#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from treeFunctions import *
import Styles
#Styles.tdrStyle()
import numpy
import ROOT
import argparse
import ConfigParser

def deltaPhi( phi1, phi2):
	import math
	result = phi1 - phi2
	while result > math.pi:
		result -= 2*math.pi
	while result <= -math.pi:
		result += 2*math.pi
	return result

def deltaR( object1, object2 ):
	import math
	return math.sqrt( (object1.eta-object2.eta)**2 + (deltaPhi(object1.phi,object2.phi))**2 )

def gammaSelectionClone( tree, newTreeName ):
	""" Clones a tree and get access to photon vector"""
	newTree = tree.CloneTree(0)
	newTree.SetName( newTreeName )
	photons = ROOT.std.vector("tree::Photon")()
	newTree.SetBranchAddress("photon", photons )

	return newTree, photons

def combineVectors( v1, v2 ):
	for element in v1:
		v2.push_back(element)
	return v2

def recElectronMatching( genElectrons, photonElectrons ):
	for recE in photonElectrons:
		minDeltaR = 50
		minPt = 8000
		minIndex = -1
		for i in range( genElectrons.size()):
			dr = deltaR( recE, genElectrons[i] )
			pt = ( recE.pt - genElectrons[i].pt ) / genElectrons[i].pt
			if dr < minDeltaR:
				minIndex = i
				minDeltaR = dr
				minPt = pt

		if minDeltaR < 0.02 and abs(minPt) < 0.5:
			photonElectrons[minIndex].phi = 5
			genElectrons[minIndex].phi = 5

	return genElectrons, photonElectrons

def genElectronMatching( photons, photonElectrons, genElectrons, hist ):
	for genE in genElectrons:
		minDeltaR = 50
		minPt = 8000
		minIndex = -1
		typeName = ""
		for i in range( photons.size()):
			dr = deltaR( genE, photons[i] )
			pt = ( photons[i].pt - genE.pt ) / genE.pt
			if dr < minDeltaR:
				minIndex = i
				minDeltaR = dr
				minPt = pt
				typeName = "g"
		for i in range( photonElectrons.size()):
			dr = deltaR( genE, photonElectrons[i] )
			pt = ( photonElectrons[i].pt - genE.pt ) / genE.pt
			if dr < minDeltaR:
				minIndex = i
				minDeltaR = dr
				minPt = pt
				typeName = "e"

		#if minDeltaR < 0.02 and abs(minPt) < 0.5:
		if minDeltaR<50:
			hist.Fill( minPt, minDeltaR )
			if typeName == "e":
				photonElectrons[minIndex].pixelseed = -2
			if typeName == "g":
				photons[minIndex].pixelseed = -2
	return photons, photonElectrons, hist

def weightTree( photonTree, photonJetTree ):
	# a tree can't be modified, therefore a new tree is saved within the same file

	# this is the binning with witch the fo reweighting will be done and has to be studied
	nBins = 10
	xMin = 80
	xMax = 800
	h_photonPt = createHistoFromTree( photonTree, "photon[0].pt", "", nBins, xMin, xMax )
	h_jetPt = createHistoFromTree( photonJetTree, "photon[0].pt", "", nBins, xMin, xMax )

	# h_jetPt is now histogram with w^{-1} = h_jetPt / h_photonPt
	h_jetPt.Divide( h_photonPt )
	h_jetPt.SetTitle(";p_{T,#gamma};w^{-1}")
	#h_jetPt.Draw()
	#raw_input()

	weightTree = ROOT.TTree("QCDWeightTree", "my TFriend for weights")
	import numpy

	# a python float corresponds to a root double
	weight = numpy.zeros(1, dtype=float)
	weightError = numpy.zeros(1, dtype=float)
	weightTree.Branch( "weight", weight, "weight/D" )
	weightTree.Branch( "weightError", weightError, "weightError/D" )

	for event in photonJetTree:
		try:
			pt = event.photon.at(0).pt
		except:
			pt = 0
		bin = h_jetPt.FindBin( pt )
		weight[0] = h_jetPt.GetBinContent( bin )
		weightError[0] = h_jetPt.GetBinError( bin )
		weightTree.Fill()

	return weightTree

def splitCandidates( inputFileName, nExpected, isQCD, isEWK):
	outputFileName = "slim"+inputFileName

	eventHisto = readHisto( inputFileName )
	nGenerated = eventHisto.GetBinContent(1)

	tree = readTree( inputFileName )

	photonTree, photons = gammaSelectionClone( tree, "photonTree" )
	photonJetTree, photonJets = gammaSelectionClone( tree, "photonJetTree" )
	photonElectronTree, photonElectrons = gammaSelectionClone( tree, "photonElectronTree" )
	#genElectronTree = gammaSelectionClone( tree, "genElectronTree" )[0]

	weight = numpy.zeros(1, dtype=float)
	photonTree.SetBranchAddress("weight", weight)
	photonElectronTree.SetBranchAddress("weight", weight)
	photonJetTree.SetBranchAddress("weight", weight)
	#genElectronTree.SetBranchAddress("weight", weight)

	#genElectrons = ROOT.std.vector("tree::Particle")()
	#photonElectronTree.SetBranchAddress("genElectron", genElectrons )


	hist2 = ROOT.TH2F( inputFileName, "", 100, -.30, .30, 100, 0, 0.06 )
	hist2.SetTitle(";#frac{p_{T,rec}-p_{T,gen}}{p_{T,gen}};min #DeltaR")

	for event in tree:
		if not event.GetReadEntry()%100000:
			print 100*event.GetReadEntry()/event.GetEntries(), "%"
		weight[0] = event.weight * nExpected / nGenerated
		photons.clear()
		photonElectrons.clear()
		photonJets.clear()

		for gamma in event.photon:

			# cuts for every object to rject spikes
			if gamma.r9 < 1 \
			and gamma.sigmaIetaIeta > 0.001 \
			and abs(gamma.eta) < 1.479 \
			and gamma.pt > 20:

				# look for gamma and electrons
				if gamma.hadTowOverEm < 0.05 \
				and gamma.sigmaIetaIeta < 0.012 \
				and gamma.chargedIso < 2.6 \
				and gamma.neutralIso < 3.5 + 0.04*gamma.pt \
				and gamma.photonIso < 1.3 + 0.005*gamma.pt:
					if gamma.pixelseed:
						photonElectrons.push_back( gamma )
					else:
						photons.push_back(gamma)

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

		photons, photonElectrons, hist2 = genElectronMatching( photons, photonElectrons, event.genElectron, hist2 )

		#genElectrons = recElectronMatching( event.genElectron, photonElectrons )

		if photons.size() > 0:
			photonTree.Fill()
		else:
			if photonElectrons.size() > 0:
				photonElectronTree.Fill()
			if photonJets.size() > 0:
				photonJetTree.Fill()
		#if (photons.size()>0 or photonElectrons.size()>0) and event.genElectron.size() > 0:
		#	genElectronTree.Fill()

	can = ROOT.TCanvas()
	can.cd()
	hist2.Draw("colz")
	can.SaveAs("pt_r_%s_new.pdf"%inName[0:5])

	if isQCD:
		qcdWeights = weightTree( photonTree, photonJetTree )


	f = ROOT.TFile( outputFileName, "recreate" )
	f.cd()
	photonTree.Write()
	photonJetTree.Write()
	photonElectronTree.Write()
	#genElectronTree.Write()
	if isQCD:
		qcdWeights.Write()
	f.Write()
	f.Close()


if __name__ == "__main__":
	# include knowledge about objects saved in the tree
	ROOT.gSystem.Load("libTreeObjects.so")

	arguments = argparse.ArgumentParser( description="Slim tree" )
	arguments.add_argument("--input", default=["WJets_V01.03_tree.root"], nargs="+" )
	opts = arguments.parse_args()

	ROOT.gROOT.SetBatch()


	datasetConf = ConfigParser.SafeConfigParser()
	datasetConf.read( "dataset.cfg" )

	integratedLumi = 19300 #pb

	for inName in opts.input:
		shortName = None
		for configName in datasetConf.sections():
			if inName.count( configName ):
				shortName = configName
				crosssection = datasetConf.getfloat( configName, "crosssection" )

		qcd = False
		ewk = False
		if shortName in ["QCD_250_500", "QCD_500_1000", "QCD_1000_inf","GJets"]:
			qcd = True
		if shortName in ["TTJets", "WJets"]:
			ewk = True

		# N = L * sigma
		nExpected = integratedLumi * crosssection

		splitCandidates( inName, nExpected, isQCD=qcd, isEWK=ewk )

