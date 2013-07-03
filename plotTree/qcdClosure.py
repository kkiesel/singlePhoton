#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
import argparse
import ConfigParser
from math import sqrt
from multiplot import *
from treeFunctions import *
import ratios
import Styles

Styles.tdrStyle()
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libTreeObjects.so")

# to use user defined help message, sys.arv has to be sent to python and not
# to TApplication
ROOT.PyConfig.IgnoreCommandLineOptions = True

def writeWeights( fileName, tree, h_weight, variable ):
	f = ROOT.TFile( fileName, "update" )

	weightTree = ROOT.TTree("weightTree", "Tree containing QCD weights" )
	import numpy
	weight = numpy.zeros( 1, dtype=float)
	weight_error = numpy.zeros( 1, dtype=float)

	weightTree.Branch( "w_qcd", weight, "w_qcd/D" )
	weightTree.Branch( "w_qcd_error", weight_error, "w_qcd_error/D" )

	variable = variable.replace("[0]",".at(0)")
	variable = "tree."+variable

	for event in tree:
		if not event.GetReadEntry()%10000:
			print "%s / %s"%(event.GetReadEntry(), event.GetEntries() )
		pt = eval( variable )
		weight[0] = h_weight.GetBinContent( h_weight.FindBin( pt ) )
		weight_error[0] = h_weight.GetBinError( h_weight.FindBin( pt ) )
		weightTree.Fill()

	f.cd()
	weightTree.Write()
	f.Close()

def clearTreeFromJets( inTree, jets, _deltaR=.3 ):
	from splitCandidatesQCDtest import deltaR
	newTree = inTree.CloneTree(0)
	newTree.SetBranchAddress( "jet", jets )
	for event in inTree:
		jets.clear()
		for j in event.jet:
			if deltaR( j, event.photon[0] ) > _deltaR:
				jets.push_back( j )
		if jets.size() > 1:
			newTree.Fill()

	# Set the size of the tree in memory to 10MB or get memory overflow.
	#newTree.SetMaxVirtualSize(int(1e7))
	return newTree

def qcdClosure( fileName, opts, cuts ):
	# process cut string for readability
	cutSaveAffix = "%s_%s_%s_%s_%s"%(cuts["chIso"], cuts["nIsoRel"], cuts["nIso"], cuts["pIsoRel"], cuts["pIso"])
	cutText = ROOT.TPaveText(0,.95, 1, .99, "ndc" )
	cutText.AddText( "#sigma_{i#etai#eta}<0.014  Iso^{#pm}<%s  Iso^{0}<%s+%sp_{T}  Iso^{#gamma}<%s+%sp_{T}"%(cuts["chIso"],cuts["nIso"],cuts["nIsoRel"],cuts["pIso"],cuts["pIsoRel"]) )
	cutText.SetFillColor(0)
	cutText.SetLineWidth(0)
	cutText.SetLineColor(0)
	cutText.SetShadowColor(0)
	import re
	datasetAffix = re.match(".*slim(.*)_V.*", fileName ).groups()[0]

	# Apply the cuts and divide in control and signal region
	signalCut = "met > 100"
	controlCut = "met <= 100"
	foCut = "Min$(photon.pt) > 80" \
			+"&& Max$(photon.sigmaIetaIeta) < 0.014" \
			+"&& Max$(photon.hadTowOverEm) < 0.05" \
			+"&& Max$(photon.chargedIso) < 15" \
			+"&& Max$(photon.neutralIso-0.04*photon.pt)<3.5" \
			+"&& Max$(photon.photonIso-0.005*photon.pt) < 1.3" \
			+"&& ( Min$(photon.chargedIso)>2.6 || Min$(photon.sigmaIetaIeta)>0.012 )" \
			+"&& @photon.size()>0"

	foCut = "Max$(photon.hadTowOverEm) < 0.05" \
			+"&& Max$(photon.chargedIso) < %s"%cuts["chIso"] \
			+"&& Max$(photon.neutralIso- %s *photon.pt)< %s "%(cuts["nIsoRel"], cuts["nIso"]) \
			+"&& Max$(photon.photonIso- %s *photon.pt) < %s "%(cuts["pIsoRel"], cuts["pIso"]) \
			+"&& @photon.size()>0" # +"&& ( Min$(photon.chargedIso)>2.6 || Min$(photon.sigmaIetaIeta)>0.012)"

	gSignalTree = readTree( fileName, "photonTree").CopyTree( signalCut )
	gControlTree = readTree( fileName, "photonTree").CopyTree( controlCut )
	foSignalTree = readTree( fileName, "photonJetTree").CopyTree( signalCut+"&&"+foCut )
	foControlTree = readTree( fileName, "photonJetTree").CopyTree( controlCut+"&&"+foCut )

	# Clean jets from fake object
	jets = ROOT.std.vector("tree::Jet")()
	foControlTree = clearTreeFromJets( foControlTree, jets )
	foSignalTree = clearTreeFromJets( foSignalTree, jets )

	can = ROOT.TCanvas()
	can.cd()

	weightVariable = "photon[0].ptJet"
	weight_numerator = getHisto( gControlTree, weightVariable, color=1 )
	weight_denominator = getHisto( foControlTree, weightVariable, color=46 )

	multiToWeight = Multihisto()
	multiToWeight.leg.SetHeader( datasetAffix )
	multiToWeight.addHisto( weight_numerator, "#gamma" )
	multiToWeight.addHisto( weight_denominator, "#gamma_{jet}" )
	multiToWeight.Draw()
	cutText.Draw()
	SaveAs(can, "plots", "qcd_pt_preWeighting_%s_%s"%(datasetAffix,cutSaveAffix) )

	weight = divideHistos( weight_numerator, weight_denominator )
	weight.GetYaxis().SetTitle("w")
	weight.Draw()
	cutText.Draw()
	SaveAs( can, "plots", "qcd_weight_%s_%s"%(datasetAffix,cutSaveAffix) )

	writeWeights( fileName, foControlTree, weight, weightVariable )

	for tree in [ foSignalTree, foControlTree ]:
		tree.AddFriend( "weightTree", fileName )

	plots = [ "photon[0].r9", "photon[0].sigmaIetaIeta", "photon[0].hadTowOverEm",
			"photon[0].chargedIso", "photon[0].neutralIso", "photon[0].photonIso",
			"electron.pt", "muon.pt", "met", "ht", "nVertex",
			"Length$(electron.pt)", "Length$(muon.pt)",
			"photon[0].eta", "photon[0].ptJet", "photon[0].pt", "Length$(photon.pt)",
			"jet[0].pt", "jet.pt", "jet[0].eta", "jet.eta", "Length$(jet.pt)",
			"jet.bCSV", "jet[0].bCSV",
			"genPhoton.pt", "Length$(genPhoton.pt)" ]

	plots = [ "met", "ht", "photon[0].ptJet", "Length$(jet.pt)" ]

	for plot in plots:
		# The first attempt to get the histogram is only to get the minimal
		# and maximal value on the x-axis, for not predefined binning
		plotFo = plot.replace( "jet", "cleanedJet" )
		h_gamma = getHisto( gSignalTree, plot )
		h_fo = getHisto( foSignalTree, plotFo, cut="@cleanedJet.size()>1" )
		xMin, xMax = getXMinXMax( [ h_gamma, h_fo ] )

		# for integers, adjust nBins and shift by 0.5
		if "Length$(" in plot or "nVertex" == plot:
			xMin -= .5
			xMax += .5
			nBins = int(xMax-xMin)
		else:
			nBins = 20

		# there is a muon with 7000 GeV, which destroys the automatic binning
		if plot == "muon.pt":
			xMax = 200

		h_gamma = getHisto( gSignalTree, plot, color=1, nBins=nBins, firstBin=xMin,lastBin=xMax )
		h_fo = getHisto( foSignalTree, plotFo, cut="@cleanedJet.size()>1", weight="weight*w_qcd", color=46, nBins=nBins, firstBin=xMin,lastBin=xMax )
		h_fo_error = getQCDErrorHisto( foSignalTree, plotFo, cut="@cleanedJet.size()>1",nBins=nBins, firstBin=xMin, lastBin=xMax )
		h_fo_error.SetFillColor(2)
		h_fo_error.SetFillStyle(3254)
		h_fo_error.SetMarkerSize(0)

		muhisto = Multihisto()
		muhisto.leg.SetHeader( datasetAffix )
		muhisto.addHisto( h_gamma, "#gamma", draw="hist e0" )
		muhisto.addHisto( h_fo, "#gamma_{jet}#upointw", draw="hist e0")
		muhisto.addHisto( h_fo_error, "#sigma_{w}", draw="e2")

		hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
		hPad.cd()
		muhisto.Draw()
		cutText.Draw()

		ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
		ratioPad.cd()
		ratioPad.SetLogy(0)
		ratioGraph = ratios.RatioGraph(h_fo, h_gamma)
		ratioGraph.draw(ROOT.gPad, yMin=0.5, yMax=1.5, adaptiveBinning=False, errors="yx")
		ratioGraph.graph.Draw("same p e0") # draw nice points
		ratioGraph.hAxis.SetYTitle( "#gamma_{pred}/#gamma")
		can.cd()
		hPad.Draw()
		ratioPad.Draw()
		SaveAs(can, "plots","qcd_%s_%s_afterWeighting_%s"%(datasetAffix, plot,cutSaveAffix))
		ROOT.SetOwnership( hPad, False )
		ROOT.SetOwnership( ratioPad, False )
	for tree in [foSignalTree, foControlTree, gSignalTree, gControlTree ]:
		tree.Delete()
		del tree


def drange( start, stop, step ):
	r = start
	while r < stop:
		yield r
		r += step

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Calculate weighting "
			+"factors for QCD background estimation." )
	arguments.add_argument("--input", default=[""], nargs="+" )
	arguments.add_argument("--plot", default="met" )
	arguments.add_argument("--save", action="store_true", help="Save canvas as pdf.")
	opts = arguments.parse_args()

	for inName in opts.input:
		cuts = { "chIso":  4.6,
				"nIsoRel": 0.05,
				"nIso":    4.5,
				"pIsoRel": 0.005,
				"pIso":    1.3 }

		qcdClosure( inName, opts, cuts )
		for nIsoRel in drange( .005, 0.024, 0.001 ):
			cuts["pIsoRel"] = nIsoRel
			for nIso in drange( 1.3, 13.5, .5 ):
				cuts["pIso"] = nIso
