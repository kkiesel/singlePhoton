#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
ROOT.gSystem.Load("libTreeObjects.so")
import argparse
from multiplot import *
from treeFunctions import *
import Styles
style = Styles.tdrStyle()
style.SetOptLogy(0)
ROOT.TGaxis.SetMaxDigits(3)
import os
from math import sqrt

import ConfigParser
axisConf = ConfigParser.SafeConfigParser()
axisConf.read("axis.cfg")

def extractHisto( dataset, plot ):
	yutarosBinning = [ 25, 35, 40, 50, 60, 80, 100 ]
	return createHistoFromTree( dataset.tree, plot, "weight*(%s)"%(dataset.additionalCut), nBins=yutarosBinning)

def yutarosHistogramMC():
	import array
	yutarosBinning = [ 25, 35, 40, 50, 60, 80, 100 ]
	yutaro = ROOT.TH1F("yutaro", "Yutaros fake rate", len(yutarosBinning)-1, array.array('d', yutarosBinning) )
	fakeRate = [(1.31, .04, .06),
			(1.46, .02, .06),
			(1.48, .01, .01),
			(1.11, .02, .02),
			(1.11, .04, .0 ),
			(0.84, .05, .02 )] # in %
	for i, val in enumerate(fakeRate):
		yutaro.SetBinContent( i+1, val[0]/100 )
		yutaro.SetBinError  ( i+1, sqrt(val[1]**2 + val[2]**2)/100)
	yutaro.SetLineColor(2)
	yutaro.SetMarkerColor(2)
	return yutaro

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="photon.pt" )
	arguments.add_argument( "--input", default="EWK_V01.12_tree.root" )
	arguments.add_argument( "--savePrefix", default="new" )
	opts = arguments.parse_args()

	ROOT.gROOT.SetBatch()
	import re
	# dataset name is from beginning till first '_'
	slimFileName = opts.input.replace( os.path.basename(opts.input), "slim"+os.path.basename(opts.input))
	datasetAffix = re.match("slim([^_]*)_.*", slimFileName ).groups()[0]

	h_gamma = extractHisto( Dataset( slimFileName, "photonTree", "photon.genInformation==1" ), opts.plot )
	h_e = extractHisto( Dataset( slimFileName, "photonElectronTree", "photon.genInformation==1" ), opts.plot )

	yuFakeRate = divideHistos( h_gamma, addHistos( [h_gamma, h_e] ) )

	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(0)

	datasetLabel = ROOT.TPaveText(.4,.94,.6,1, "ndc")
	datasetLabel.SetFillColor(0)
	datasetLabel.SetBorderSize(0)
	datasetLabel.AddText( datasetAffix )

	label, unit = readAxisConf( opts.plot, axisConf )
	yuFakeRate.SetTitle(";%s%s;#gamma / (e+#gamma)"%(label,unit) )
	yuOrig = yutarosHistogramMC()

	mhisto = Multihisto()
	mhisto.legendOption = "lp"
	mhisto.addHisto( yuFakeRate, "MC info", draw="")
	mhisto.addHisto( yuOrig, "DY tag&probe", draw="")
	mhisto.leg.SetHeader("Object matching")
	mhisto.Draw()

	datasetLabel.Draw()

	saveName = "%s_%s_%s_%s"%(yuFakeRate.GetYaxis().GetTitle(),datasetAffix,opts.plot,opts.savePrefix)
	saveName = saveName.replace("/","VS")
	saveName = saveName.replace(" ","_")
	unallowedCharacters = ["{","}","(",")","#","|","."]
	for char in unallowedCharacters:
		saveName = saveName.replace( char, "" )

	can.SaveAs("plots/%s.pdf"%saveName)


