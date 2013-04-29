#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from treeFunctions import *
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

def genMatching( photons, photonElectrons, genElectrons, hist ):
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

		if minIndex != -1:
			hist.Fill( minPt, minDeltaR )
			if typeName == "e":
				photonElectrons[minIndex].pixelseed = -2
			if typeName == "g":
				photons[minIndex].pixelseed = -2
	return photons, photonElectrons, hist

def process( inputFileName, nExpected ):
	outputFileName = "slim"+inputFileName

	eventHisto = readHisto( inputFileName )
	nGenerated = eventHisto.GetBinContent(1)

	tree = readTree( inputFileName )

	photonTree, photons = gammaSelectionClone( tree, "photonTree" )
	photonJetTree, photonJets = gammaSelectionClone( tree, "photonJetTree" )
	photonElectronTree, photonElectrons = gammaSelectionClone( tree, "photonElectronTree" )

	pu_weight = numpy.zeros(1, dtype=float)
	ROOT.gROOT.ProcessLine("float pu_weight;")
	from ROOT import pu_weight

	print type( pu_weight )
	photonTree.SetBranchAddress("pu_weight", pu_weight)

	hist2 = ROOT.TH2F("blub2", "", 100, -.2, 1, 100, 0, 0.06 )
	hist2.SetTitle(";#frac{p_{T,rec}-p_{T,gen}}{p_{T,gen}};min #DeltaR")

	for event in tree:
		pu_weight[0] = event.pu_weight * nExpected / nGenerated
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
				if gamma.pt > 80 \
				and gamma.hadTowOverEm < 0.05 \
				and gamma.sigmaIetaIeta < 0.014 \
				and gamma.chargedIso < 15 \
				and gamma.neutralIso < 3.5 + 0.04*gamma.pt \
				and gamma.photonIso < 1.3 + 0.005*gamma.pt \
				and not gamma.pixelseed \
				and ( gamma.sigmaIetaIeta>=0.012 or gamma.chargedIso>=2.6):
					photonJets.push_back( gamma )

		photons, photonElectrons, hist2 = genMatching( photons, photonElectrons, event.genElectron, hist2 )


		if photons.size() > 0:
			photonTree.Fill()
		else:
			if photonElectrons.size() > 0:
				photonElectronTree.Fill()
			if photonJets.size() > 0:
				photonJetTree.Fill()

	#hist2.Draw("colz")
	#raw_input()

	eventHisto = readHisto( inputFileName )

	f = ROOT.TFile( outputFileName, "recreate" )
	f.cd()
	photonTree.Write()
	photonJetTree.Write()
	photonElectronTree.Write()
	f.Write()
	f.Close()


if __name__ == "__main__":
	# include knowledge about objects saved in the tree
	ROOT.gSystem.Load("libTreeObjects.so")

	arguments = argparse.ArgumentParser( description="Slim tree" )
	arguments.add_argument("--input", default=["WJets_V01.03_tree.root"], nargs="+" )
	opts = arguments.parse_args()



	datasetConf = ConfigParser.SafeConfigParser()
	datasetConf.read( "dataset.cfg" )

	integratedLumi = 19300 #pb

	for inName in opts.input:
		for configName in datasetConf.sections():
			if inName.count( configName ):
				crosssection = datasetConf.getfloat( configName, "crosssection" )

		# N = L * sigma
		nExpected = integratedLumi * crosssection

		process( inName, nExpected )

