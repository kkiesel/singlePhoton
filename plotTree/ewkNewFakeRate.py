#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
import re
import os
import argparse
from multiplot import *
from treeFunctions import *
from math import sqrt

import Styles
style = Styles.tdrStyle()
style.SetOptLogy(0)
ROOT.TGaxis.SetMaxDigits(3)
ROOT.gROOT.SetBatch()

import ConfigParser
axisConf = ConfigParser.SafeConfigParser()
axisConf.read("axis.cfg")

ROOT.gSystem.Load("libTreeObjects.so")

def extractHisto( dataset, plot ):
	yutarosBinning = [ 25, 35, 40, 50, 60, 80, 100 ]
	return createHistoFromTree( dataset.tree, plot, "weight*(%s)"%(dataset.additionalCut), nBins=yutarosBinning)

def yutarosHistogramMC( color=2 ):
	"""Creates Yutaros fake rate histogram for mc"""
	import array
	yutarosBinning = [ 25, 35, 40, 50, 60, 80, 100 ]
	yutaro = ROOT.TH1F("yutaro", "Yutaros fake rate", len(yutarosBinning)-1, array.array('d', yutarosBinning) )
	# values copied from presentation
	fakeRate = [(1.31, .04, .06),
			(1.46, .02, .06),
			(1.48, .01, .01),
			(1.11, .02, .02),
			(1.11, .04, .0 ),
			(0.84, .05, .02 )]
	# set bin values and errors
	for i, val in enumerate(fakeRate):
		yutaro.SetBinContent( i+1, val[0]/100 )
		yutaro.SetBinError  ( i+1, sqrt(val[1]**2 + val[2]**2)/100)
	yutaro.SetLineColor( color )
	yutaro.SetMarkerColor( color )
	return yutaro

def drawDatasetLabel( datasetAffix ):
	"""Draws sample info on top of the canvas."""
	datasetLabel = ROOT.TPaveText(.4,.94,.6,1, "ndc")
	datasetLabel.SetFillColor(0)
	datasetLabel.SetBorderSize(0)
	datasetLabel.AddText( datasetAffix )
	return datasetLabel

def manipulateSaveName( saveName ):
	"""Replace some charakters, so root nor unix have problems to read them."""
	saveName = saveName.replace("/","VS")
	saveName = saveName.replace(" ","_")
	unallowedCharacters = ["{","}","(",")","#","|","."]
	for char in unallowedCharacters:
		saveName = saveName.replace( char, "" )
	return saveName

def plotNewFakeRate( fileName, opts ):
	# dataset name is from beginning till first '_'
	slimFileName = fileName.replace( os.path.basename(fileName), "slim"+os.path.basename(fileName))
	datasetAffix = re.match("slim([^_]*)_.*", slimFileName ).groups()[0]

	h_gamma = extractHisto( Dataset( slimFileName, "photonTree", "photon.isGenElectron()" ), opts.plot )
	h_e = extractHisto( Dataset( slimFileName, "photonElectronTree", "photon.isGenElectron()" ), opts.plot )

	fakeRate = divideHistos( h_gamma, addHistos( [h_gamma, h_e] ) )

	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(0)

	label, unit, binning = readAxisConf( opts.plot, axisConf )
	fakeRate.SetTitle(";%s%s;#gamma / (e+#gamma)"%(label,unit) )
	yuOrig = yutarosHistogramMC()

	mhisto = Multihisto()
	mhisto.legendOption = "lp"
	mhisto.addHisto( fakeRate, "MC info", draw="")
	mhisto.addHisto( yuOrig, "DY tag&probe", draw="")
	mhisto.leg.SetHeader("Object matching")
	mhisto.Draw()

	label = drawDatasetLabel( datasetAffix )
	label.Draw()

	saveName = "%s_%s_%s_%s"%(fakeRate.GetYaxis().GetTitle(),datasetAffix,opts.plot,opts.savePrefix)
	can.SaveAs( "plots/%s.pdf"%manipulateSaveName(saveName) )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="photon.pt" )
	arguments.add_argument( "--input", default=["WJets_V01.12_tree.root"], nargs="+" )
	arguments.add_argument( "--savePrefix", default="new" )
	opts = arguments.parse_args()

	for inName in opts.input:
		plotNewFakeRate( inName, opts )

