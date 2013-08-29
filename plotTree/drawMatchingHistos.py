#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
import argparse
import ConfigParser
import re
from math import sqrt
from multiplot import *
from treeFunctions import *
import ratios
import Styles

Styles.tdrStyle2D()
ROOT.gSystem.Load("libTreeObjects.so")

# to use user defined help message, sys.arv has to be sent to python and not
# to TApplication
ROOT.PyConfig.IgnoreCommandLineOptions = True

def makePlot( inName, histoName, drawCuts ):
	datasetAffix = re.match("(slim)?(.*)_V.*", inName ).groups()[-1]
	h = readHisto( inName, histoName )

	can = ROOT.TCanvas()
	can.cd()
	h.Draw("colz")

	if drawCuts:
		drCut = .3
		relPtCut = .95

		drCutLine = ROOT.TLine( drCut, relPtCut, drCut, 4 )
		relPtCutLine = ROOT.TLine( 0, relPtCut, drCut, relPtCut )
		for line in drCutLine, relPtCutLine:
			line.SetLineWidth(2)
			line.SetLineColor(2)
			line.Draw()

	can.SaveAs( "plots/matchingPhotonJet_%s_%s.pdf"%(datasetAffix, histoName) )
	return h

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Calculate weighting "
			+"factors for QCD background estimation." )
	arguments.add_argument("--input", default=[], nargs="+" )
	opts = arguments.parse_args()

	for inName in opts.input:
		drawCuts = False
		f = ROOT.TFile( inName )
		tList = f.GetListOfKeys()
		for obj in tList:
			histoName = obj.GetName()
			if  "matchingPhotonJet" in histoName:
				drawCuts = True

			if "matching" in histoName:
				makePlot( inName, histoName, drawCuts )

	try:
		h = []
		for inName in opts.input:
			h.append( makePlot( inName, "matchingPhotonJetdRPtGamma", True ) )
			h.append( makePlot( inName, "matchingPhotonJetdRPt", True ) )
		allHistos = addHistos( h )
		can = ROOT.TCanvas()
		can.cd()
		allHistos.Draw("colz")

		drCut = .3
		relPtCut = .95

		drCutLine = ROOT.TLine( drCut, relPtCut, drCut, 4 )
		relPtCutLine = ROOT.TLine( 0, relPtCut, drCut, relPtCut )
		for line in drCutLine, relPtCutLine:
			line.SetLineWidth(2)
			line.SetLineColor(2)
			line.Draw()

		can.SaveAs( "plots/matchingPhotonJet.pdf" )
	except:
		pass

