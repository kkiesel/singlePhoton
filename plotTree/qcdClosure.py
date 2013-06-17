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

def writeWeights( fileName, tree, h_weight ):
	f = ROOT.TFile( fileName, "update" )

	weightTree = ROOT.TTree("weightTree", "Tree containing QCD weights" )
	import numpy
	weight = numpy.zeros( 1, dtype=float)
	weight_error = numpy.zeros( 1, dtype=float)

	weightTree.Branch( "w_qcd", weight, "w_qcd/D" )
	weightTree.Branch( "w_qcd_error", weight_error, "w_qcd_error/D" )

	for event in tree:
		if not event.GetReadEntry()%10000:
			print "%s / %s"%(event.GetReadEntry(), event.GetEntries() )
		pt = tree.photon.at(0).ptJet
		weight[0] = h_weight.GetBinContent( h_weight.FindBin( pt ) )
		weight_error[0] = h_weight.GetBinError( h_weight.FindBin( pt ) )
		weightTree.Fill()

	f.cd()
	weightTree.Write()
	f.Close()

def qcdClosure( fileName, opts ):
	signalCut = "met > 100"
	controlCut = "met <= 100"
	foCut = "Min$(photon.pt) > 80" \
			+"&& Max$(photon.sigmaIetaIeta) < 0.014" \
			+"&& Max$(photon.hadTowOverEm) < 0.05" \
			+"&& Max$(photon.chargedIso) < 15" \
			+"&& Max$(photon.neutralIso-0.04*photon.pt)<3.5" \
			+"&& Max$(photon.photonIso-0.005*photon.pt) < 1.3" \
			+"&& ( Min$(photon.chargedIso)>2.6 || Min$(photon.sigmaIetaIeta)>0.012) " \
			+"&& @photon.size()>0"

	gSignalTree = readTree( fileName, "photonTree").CopyTree( signalCut )
	gControlTree = readTree( fileName, "photonTree").CopyTree( controlCut )
	foSignalTree = readTree( fileName, "photonJetTree").CopyTree( signalCut+"&&"+foCut )
	foControlTree = readTree( fileName, "photonJetTree").CopyTree( controlCut+"&&"+foCut )

	can = ROOT.TCanvas()
	can.cd()

	weight_numerator = getHisto( gControlTree, "photon[0].ptJet", color=1 )
	weight_denominator = getHisto( foControlTree, "photon[0].ptJet", color=6 )

	multiToWeight = Multihisto()
	multiToWeight.addHisto( weight_numerator, "#gamma" )
	multiToWeight.addHisto( weight_denominator, "#gamma_{jet}" )
	multiToWeight.Draw()
	SaveAs(can, "plots", "qcd_pt_preWeighting")

	weight = divideHistos( weight_numerator, weight_denominator )
	weight.GetYaxis().SetTitle("w")
	weight.Draw()
	SaveAs( can, "plots", "qcd_weight" )

	writeWeights( fileName, foControlTree, weight )

	for tree in [ foSignalTree, foControlTree ]:
		tree.AddFriend( "weightTree", fileName )

	for plot in ["photon[0].ptJet", "met", "ht", "jet[0].pt", "nVertex"]:
		# The first attempt to get the histogram is only to get the minimal
		# and maximal value on the x-axis, for not predefined binning
		h_gamma = getHisto( gControlTree, plot )
		h_fo = getHisto( foControlTree, plot )
		xMin, xMax = getXMinXMax( [ h_gamma, h_fo ] )

		h_gamma = getHisto( gControlTree, plot, color=1, firstBin=xMin,lastBin=xMax )
		h_fo = getHisto( foControlTree, plot, weight="weight*w_qcd", color=6,firstBin=xMin,lastBin=xMax )
		h_fo_error = getQCDErrorHisto( foControlTree, plot, firstBin=xMin, lastBin=xMax )
		h_fo_error.SetFillColor(2)
		h_fo_error.SetFillStyle(3254)
		h_fo_error.SetMarkerSize(0)

		muhisto = Multihisto()
		muhisto.addHisto( h_gamma, "#gamma", draw="hist e0" )
		muhisto.addHisto( h_fo, "#gamma_{jet}#upointw", draw="hist e0")
		muhisto.addHisto( h_fo_error, "#sigma_{w}", draw="e2")

		hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
		hPad.cd()
		muhisto.Draw()

		ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
		ratioPad.cd()
		ratioPad.SetLogy(0)
		ratioGraph = ratios.RatioGraph(h_fo, h_gamma)
		ratioGraph.draw(ROOT.gPad, yMin=0.5, yMax=1.5, adaptiveBinning=False, errors="yx")
		ratioGraph.hAxis.SetYTitle( "#gamma_{pred}/#gamma")
		can.cd()
		hPad.Draw()
		ratioPad.Draw()
		SaveAs(can, "plots","qcd_%s_afterWeighting"%plot )
		ROOT.SetOwnership( hPad, False )
		ROOT.SetOwnership( ratioPad, False )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Calculate weighting "
			+"factors for QCD background estimation." )
	arguments.add_argument("--input", default=[""], nargs="+" )
	arguments.add_argument("--plot", default="met" )
	arguments.add_argument("--save", action="store_true", help="Save canvas as pdf.")
	opts = arguments.parse_args()

	for inName in opts.input:
		qcdClosure( inName, opts )
