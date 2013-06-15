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

def writeWeights( fileName, h_weight ):
	import array
	weight = array.array( "f", [0] )
	f = ROOT.TFile( fileName, "update" )
	tree = f.Get("photonJetTree")
	tree.Branch( "wqcd", weight, "wqcd/F" )

	for event in tree:
		if not event.GetReadEntry()%1000:
			print "%s / %s"%(event.GetReadEntry(), event.GetEntries() )
		pt = tree.photon.at(0).pt
		weight = h_weight.GetBinContent( h_weight.FindBin( pt ) )
		tree.Fill()

	tree.Write()
	f.Close()

def qcdClosure( fileName, opts ):
	commonCut = "met <= 100"
	foCut = "&& Max$(photon.sigmaIetaIeta) < 0.014" \
			+"&& Max$(photon.hadTowOverEm) < 0.05" \
			+"&& Max$(photon.chargedIso) < 15" \
			+"&& Max$(photon.neutralIso-0.0001*photon.pt)<3.5" \
			+"&& Max$(photon.photonIso-0.0001*photon.pt) < 1.3"

	gamma_set = Dataset( fileName, "photonTree", commonCut, color=1 )
	fo_set = Dataset( fileName, "photonJetTree", commonCut+foCut, color=6 )

	pt = extractHisto( gamma_set, "photon[0].ptJet" )
	pt_fo = extractHisto( fo_set, "photon[0].ptJet" )

	multihisto = Multihisto()
	multihisto.addHisto( pt, "#gamma", draw="hist" )
	multihisto.addHisto( pt_fo, "#gamma_{jet}", draw="hist" )
	can = ROOT.TCanvas()
	can.cd()
	multihisto.Draw()
	can.SaveAs("plots/qcd_pt_preWeighting.pdf")

	weight = divideHistos( pt, pt_fo )
	weight.GetYaxis().SetTitle("w")
	weight.Draw()
	can.SaveAs("plots/qcd_weight.pdf")

	writeWeights( fileName, weight )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Calculate weighting "
			+"factors for QCD background estimation." )
	arguments.add_argument("--input", default=[""], nargs="+" )
	arguments.add_argument("--plot", default="met" )
	arguments.add_argument("--save", action="store_true", help="Save canvas as pdf.")
	opts = arguments.parse_args()

	for inName in opts.input:
		qcdClosure( inName, opts )
