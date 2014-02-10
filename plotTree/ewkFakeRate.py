#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from math import sqrt

from multiplot import Multihisto
from treeFunctions import *

ROOT.TGaxis.SetMaxDigits(3)


combinedDatasets = {
		"GJets": ["GJets_200_400", "GJets_400_inf"],
		"QCD": ["QCD_250_500", "QCD_500_1000", "QCD_1000_inf"],
		"W": ["WJets_250_300", "WJets_300_400", "WJets_400_inf"],
		"WGamma": ["WGamma_50_130", "WGamma_130_inf"]
	}

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


def getFakeRateHisto( filenames, plot, cut, color ):

	gHisto, eHisto = None, None

	for fileName in filenames:
		gTree = readTree( fileName, "photonTree" )
		eTree = readTree( fileName, "photonElectronTree" )
		h_gamma = getHisto( gTree, plot, cut+"&& photons.isGen(1)", color=color )
		h_e = getHisto( eTree, plot, cut+"&& photons.isGen(1)", color=color )
		if gHisto:
			gHisto.Add( h_gamma )
			eHisto.Add( h_e )
		else:
			gHisto = h_gamma
			eHisto = h_e

	fakeRate = divideHistos( gHisto, addHistos( [gHisto, eHisto] ) )
	fakeRate.GetYaxis().SetTitle("f_{e#rightarrow#gamma}")

	return fakeRate


def plotFakeRate( filenames, plot ):

	abbrAndFile = [ (getDatasetAbbr(f), f) for f in filenames ]

	abbrMergedFiles = []

	# combine MC datasets
	sumBinned = True
	if sumBinned:
		for combiAbbr, abbrList in combinedDatasets.iteritems():
			if set(abbrList).issubset( set( [ a for a, f in abbrAndFile ] ) ):
				#print "a valid combination was found"
				abbrMergedFiles.append( (combiAbbr, [ f for a, f in abbrAndFile if a in abbrList ] )  )
				abbrAndFile = [ (a,f) for a,f in abbrAndFile if a not in abbrList ]
		for a,f in abbrAndFile:
			abbrMergedFiles.append( (a,[f]) )

	mhisto = Multihisto()
	#mhisto.setMinimum(0)
	#mhisto.setMaximum(0.02)
	mhisto.legendOption = "lp"
	#mhisto.leg.SetHeader("Object matching")

	for iColor, (abbr, filenames) in enumerate(abbrMergedFiles):
		mhisto.addHisto( getFakeRateHisto( filenames, plot, "photons.pt>35", iColor+2 ), datasetToLatex(abbr), draw="" )

	if plot == "photons.pt":
		yuOrig = yutarosHistogramMC(1)
		mhisto.addHisto( yuOrig, "DY tag&probe", draw="")

	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(0)
	mhisto.Draw()
	saveName = "%s_%s"%("fakeRate",plot)
	can.SaveAs( "plots/%s.pdf"%manipulateSaveName(saveName) )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="photons.pt" )
	opts = arguments.parse_args()

	plotFakeRate( opts.filenames, opts.plot )

