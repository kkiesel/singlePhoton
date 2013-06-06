#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
import os
import re
import argparse
from math import sqrt
from multiplot import *
from treeFunctions import *

import Styles
style = Styles.tdrStyle()
style.SetOptLogy(0)
ROOT.TGaxis.SetMaxDigits(3)
ROOT.gROOT.SetBatch()

import ConfigParser
axisConf = ConfigParser.SafeConfigParser()
axisConf.read("axis.cfg")

ROOT.gSystem.Load("libTreeObjects.so")

def yutarosHistogramMC( color=2 ):
	"""Creates Yutaros fake rate histogram for mc"""
	import array
	yutarosBinning = [ 25, 35, 40, 50, 60, 80, 100 ]
	yutaro = ROOT.TH1F("yutaro", "Yutaros fake rate", len(yutarosBinning)-1, array.array('d', yutarosBinning) )
	# fake-rate, stat and syst error in [%], copied from
	# https://twiki.cern.ch/twiki/pub/CMS/SusyPhotonID/egamma_fakerate.pdf
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

def extractHisto( dataset, plot ):
	yutarosBinning = [ 25, 35, 40, 50, 60, 80, 100 ]
	histo = createHistoFromTree( dataset.tree, plot, "weight*(%s)"%(dataset.additionalCut), nBins=yutarosBinning)
	label, unit, binning = readAxisConf( opts.plot, axisConf )
	if unit:
		unit = " ["+unit+"]"
	histo.SetTitle(";%s%s;"%(label,unit))
	return histo

def plotFakeRate( fileName, opts ):
	# dataset name is from beginning till first '_'
	slimFileName = fileName.replace( os.path.basename(fileName), "slim"+os.path.basename(fileName))
	datasetAffix = re.match("slim([^_]*)_.*", slimFileName ).groups()[0]

	genE = extractHisto( Dataset( slimFileName, "genElectronTree" ), opts.plot )
	genE_with_match = extractHisto( Dataset( slimFileName, "genElectronTree", "genElectron.phi > 4.9" ), opts.plot )
	gamma = extractHisto( Dataset( slimFileName, "photonTree" ), opts.plot )
	gamma_with_match = extractHisto( Dataset( slimFileName, "photonTree", "photon.isGenElectron()" ), opts.plot )

	label, unit, binning = readAxisConf( opts.plot, axisConf )
	e_match_e_reco = divideHistos( gamma_with_match, genE )
	e_match_e_reco.GetYaxis().SetTitle("#varepsilon_{match}#upointf_{e_{gen}#rightarrow#gamma}")

	e_match = divideHistos( genE_with_match, gamma )
	e_match.GetYaxis().SetTitle("#varepsilon_{match}")

	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(0)
	label = drawDatasetLabel( datasetAffix )

	e_match_e_reco.Draw("e")
	label.Draw()
	can.SaveAs("plots/%sEfficiencyFakeRate.pdf"%datasetAffix)

	e_match.Draw("e")
	label.Draw()
	can.SaveAs("plots/%sEfficiency.pdf"%datasetAffix)

	yutaro = yutarosHistogramMC()
	myFakeRate = divideHistos( e_match_e_reco, e_match )
	myFakeRate.GetYaxis().SetTitle("f_{e_{gen}#rightarrow #gamma}")
	myFakeRate.SetMaximum(max(myFakeRate.GetMaximum(),yutaro.GetMaximum())+0.002)
	myFakeRate.SetMinimum(min(myFakeRate.GetMinimum(),yutaro.GetMinimum())-0.002)
	myFakeRate.Draw("e")
	yutaro.Draw("same e0")

	leg = myLegend(.5,.70,.95,.92)
	leg.AddEntry( myFakeRate, myFakeRate.GetYaxis().GetTitle(), "lp")
	leg.AddEntry( yutaro, "f_{e#rightarrow#gamma} with DY", "lp" )
	leg.SetBorderSize(1)
	leg.Draw()

	label.Draw()

	saveName = "%s_%s_%s_%s"%(myFakeRate.GetYaxis().GetTitle(),datasetAffix,opts.plot,opts.savePrefix)
	can.SaveAs( "plots/%s.pdf"%manipulateSaveName(saveName) )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="photon.pt" )
	arguments.add_argument( "--input", default=["EWK_V01.12_tree.root"], nargs="+" )
	arguments.add_argument( "--savePrefix", default="new" )
	opts = arguments.parse_args()

	for inName in opts.input:
		plotFakeRate( inName, opts )


