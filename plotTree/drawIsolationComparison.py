#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

colors = {
		"GJets": ROOT.kCyan-7,
		"GJets_200_400": ROOT.kCyan-10,
		"GJets_400_inf":ROOT.kCyan-7,

		"QCD": ROOT.kCyan+3,
		"QCD_250_500": ROOT.kCyan+3,
		"QCD_500_1000": ROOT.kCyan+2,
		"QCD_1000_inf": ROOT.kCyan+1,

		"WJets": ROOT.kGreen,
		"WJets_250_inf": ROOT.kGreen,
		"WJets_250_300": ROOT.kGreen,
		"WJets_300_400": ROOT.kGreen+1,
		"WJets_400_inf": ROOT.kGreen+2,

		"TTJets": ROOT.kGreen-3,
		"TTFull": ROOT.kGreen-2,
		"TTSemi": ROOT.kGreen-1,
		"TTHadronic": ROOT.kGreen-0,
		"WGamma": ROOT.kRed-7,
		"WGamma_50_130": ROOT.kRed-7,
		"WGamma_130_inf": ROOT.kRed-6,
		"WGamma_50_inf": ROOT.kRed-6,
		"ZGamma": ROOT.kRed-4,
		"ZGammaNuNu": ROOT.kRed-3,
		"ZGammaLL": ROOT.kRed-2,
		"TTGamma": ROOT.kRed,
		}

combinedDatasets = {
		"GJets": ["GJets_200_400", "GJets_400_inf"],
		"QCD": ["QCD_250_500", "QCD_500_1000", "QCD_1000_inf"],
		"W": ["WJets_250_300", "WJets_300_400", "WJets_400_inf"],
		"WGamma": ["WGamma_50_130", "WGamma_130_inf"]
	}

def drawStackedBackground( plot, treeName, listOfFiles, sumBinned, order=False ):

	cut = "!@electrons.size() && !@muons.size()"

	# if a data histogram is present, the mc integral will scaled to data
	mcHists = []
	dataHists = []
	for fileName in listOfFiles:
		datasetAbbr = getDatasetAbbr( fileName )
		try:
			color = colors[datasetAbbr]
		except:
			color = 1
		tree = readTree( fileName, treeName )
		histo = getHisto( tree, plot, color=color, cut=cut )
		if "PhotonHad" in fileName:
			dataHists.append( (datasetAbbr, histo ) )
		else:
			mcHists.append( ( datasetAbbr, histo ) )

	# merge datahists if there
	dataHist = addHistos( [histo for datasetAbbr, histo in dataHists] ) if dataHists else None

	scale = 1.*dataHist.Integral()/addHistos( [histo for datasetAbbr, histo in mcHists] ).Integral() if dataHists else 1
	scale = 1.

	# combine MC datasets
	if sumBinned:
		for combiAbbr, abbrList in combinedDatasets.iteritems():
			if set(abbrList).issubset( set( [ a for a, h in mcHists ] ) ):
				#print "a valid combination was found"
				histosToAdd = [ h for a,h in mcHists if a in abbrList ]
				thisSum = addHistos( histosToAdd )
				mcHists = [ (a,h) for a,h in mcHists if a not in abbrList]+[(combiAbbr,thisSum)]


	mh = Multihisto()
	mh.setMinimum(0.01)
	mh.leg.SetX1NDC(0.5)
	mh.leg.SetY1NDC(0.5)
	mh.orderByIntegral = order
	if dataHist:
		mh.addHisto( dataHist, "Data", draw="e0" )

	for abbr, hist in mcHists:
		hist.Scale( scale )
		mh.addHisto( hist, datasetToLatex(abbr), toStack=True, draw="hist" )

	# add signal histos
	#signalFiles = ["slimW_1000_1020_375_V02.44_tree.root", "slimW_1200_1120_375_V02.44_tree.root"]
	signalFiles = ["slimW_1200_1120_375_V02.44_tree.root"]
	signalFiles = []
	for iColor, sf in enumerate(signalFiles):
		tree = readTree( sf, treeName )
		histo = getHisto( tree, plot, color=ROOT.kMagenta+iColor, cut=cut )
		mh.addHisto( histo, datasetToLatex(getDatasetAbbr(sf)), draw="hist" )

	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()
	info = PlotCaption(treeName=treeName)
	info.Draw()
	if dataHist:
		from myRatio import Ratio
		den = mh.stack.GetStack().Last().Clone( randomName() )
		den.SetLineColor(2)
		r = Ratio( "Data/Sim", dataHist, den )
		r.draw(0,2)

	allDatasetAbbr = getSaveNameFromDatasets( listOfFiles )
	SaveAs( can, "stackedHisto_%s_%s_%s"%(treeName, plot,allDatasetAbbr) )

def getIsolationHist( files, histName ):
	hist = None
	for f in files:
		h2 = readHisto( f, histName )
		h = h2.ProjectionY(randomName())
		if hist: hist.Add( h )
		else: hist = h
	return hist

def drawComparison( histName, mc, data ):
	hdata = getIsolationHist( data, histName )
	hdata.SetLineColor(1)


	hmc = getIsolationHist( mc, histName )
	hmc.SetLineColor(2)
	for h in hmc, hdata:
		h.SetLineWidth(2)

	hmc.Scale( hdata.Integral() / hmc.Integral() )

	mh = Multihisto()
	mh.addHisto( hdata, "Data", draw="e0" )
	mh.addHisto( hmc, "MC", draw="hist", toStack=True )
	c = ROOT.TCanvas()
	mh.Draw()

	from myRatio import Ratio
	r = Ratio( "Data/Bkg", hdata, hmc )
	r.draw(0,2)

	SavePad("isolationComparison_%s"%histName )

import re
if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "--order", action="store_true" )
	arguments.add_argument( "--sum", action="store_true" )
	opts = arguments.parse_args()

	version = 44
	g1 = "slimGJets_200_400_V02.%s_tree.root"%version
	g2 = "slimGJets_400_inf_V02.%s_tree.root"%version
	q1 = "slimQCD_250_500_V02.%s_tree.root"%version
	q2 = "slimQCD_500_1000_V02.%s_tree.root"%version
	q3 = "slimQCD_1000_inf_V02.%s_tree.root"%version

	tt = "slimTTJets_V02.%s_tree.root"%version
	w1 = "slimWJets_250_300_V02.%s_tree.root"%version
	w2 = "slimWJets_300_400_V02.%s_tree.root"%version
	w3 = "slimWJets_400_inf_V02.%s_tree.root"%version

	wg1 = "slimWGamma_50_130_V02.%s_tree.root"%44
	wg2 = "slimWGamma_130_inf_V02.%s_tree.root"%44
	tg = "slimTTGamma_V02.%s_tree.root"%44
	zgn = "slimZGammaNuNu_V02.%s_tree.root"%43
	zgl = "slimZGammaLL_V02.%s_tree.root"%43

	mc = [g1,g2,q1,q2,q3,tt,w1,w2,w3,wg1,wg2,tg,zgn,zgl]

	versionData=44
	data = [ "PhotonHad%s_V02.%s_tree.root"%(x,versionData) for x in [ "A", "B", "C", "D" ] ]



	histNames = []

	f = ROOT.TFile( g1 )
	for item in f.GetListOfKeys():
		if item.GetName() == "eventNumbers":
			continue # This histogram is used to scale
		elif re.match("TH2.", item.GetClassName() ) and "met" in item.GetName():
			histNames.append( item.GetName() )

	for histName in histNames:
		drawComparison( histName, mc, data )


