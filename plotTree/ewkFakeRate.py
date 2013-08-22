#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from math import sqrt

from multiplot import Multihisto
from treeFunctions import *

ROOT.TGaxis.SetMaxDigits(3)

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

def getFakeRateHisto( fileName, opts, color ):
	datasetAbbr = getDatasetAbbr( fileName, slim=False )

	commonCut = ""
	gTree = readTree( fileName, "photonTree" )
	eTree = readTree( fileName, "photonElectronTree" )
	h_gamma = getHisto( gTree, opts.plot, commonCut+" !photons.isGen(0)", color=color )
	h_e = getHisto( eTree, opts.plot, commonCut+" !photons.isGen(0)", color=color )

	fakeRate = divideHistos( h_gamma, addHistos( [h_gamma, h_e] ) )
	fakeRate.GetYaxis().SetTitle("f_{e#rightarrow#gamma}")

	return fakeRate


def plotFakeRate( filenames, opts ):
	mhisto = Multihisto()
	mhisto.legendOption = "lp"
	#mhisto.leg.SetHeader("Object matching")

	for iColor, filename in enumerate(filenames):
		mhisto.addHisto( getFakeRateHisto( filename, opts, iColor+2 ), getDatasetAbbr(filename,slim=False), draw="" )

	if opts.plot == "photons.pt":
		yuOrig = yutarosHistogramMC(1)
		mhisto.addHisto( yuOrig, "DY tag&probe", draw="")

	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(0)
	mhisto.Draw()
	saveName = "%s_%s_%s"%("fakeRate",opts.plot,opts.savePrefix)
	can.SaveAs( "plots/%s.pdf"%manipulateSaveName(saveName) )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="photons.pt" )
	arguments.add_argument( "--savePrefix", default="new" )
	opts = arguments.parse_args()

	plotFakeRate( opts.filenames, opts )

