#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

colors = {
		"GJets_200_400": ROOT.kCyan-10,
		"GJets_400_inf":ROOT.kCyan-7,
		"QCD_250_500": ROOT.kCyan+3,
		"QCD_500_1000": ROOT.kCyan+2,
		"QCD_1000_inf": ROOT.kCyan+1,
		"WJets": ROOT.kGreen,
		"TTJets": ROOT.kGreen-3,
		"WGamma": ROOT.kRed-7,
		"ZGamma": ROOT.kRed-4
		}


def drawStackedBackground( plot, treeName, listOfFiles, order=False ):
	# if a data histogram is present, the mc integral will scaled to data
	cut = "photons.isGen(0) && @photons.size()"
	cut = "1"
	mcHists = []
	dataHist = None
	for fileName in listOfFiles:
		datasetAbbr = getDatasetAbbr( fileName )
		try:
			color = colors[datasetAbbr]
		except:
			color = 1
		tree = readTree( fileName, treeName )
		histo = getHisto( tree, plot, color=color, cut=cut )
		if not "PhotonHad" in fileName:
			mcHists.append( (datasetAbbr, histo ) )
		else:
			dataHist = histo

	try:
		scale = 1.*dataHist.Integral()/addHistos( mcHists ).Integral()
	except:
		scale = 1.

	mh = Multihisto()
	mh.orderByIntegral = order
	if dataHist:
		mh.addHisto( dataHist, "Data", draw="e0" )

	for abbr, hist in mcHists:
		hist.Scale( scale )
		mh.addHisto( hist, abbr, toStack=True, draw="hist" )

	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()

	info = PlotCaption()
	if treeName == "photonTree":
		info.appendFront( "#gamma    " )
	elif treeName == "photonJetTree":
		info.appendFront( "#gamma_{jet}    " )
	elif treeName == "photonElectronTree":
		info.appendFront( "#gamma_{e}    " )
	else:
		info.appendFront( "%s    "%treeName )
	info.Draw()

	SaveAs( can, "stackedHisto_%s_%s_%s"%(treeName, plot,len(listOfFiles)))

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--tree", default="photonTree" )
	arguments.add_argument( "--order", action="store_true" )
	opts = arguments.parse_args()

	drawStackedBackground( opts.plot, opts.tree, opts.filenames, opts.order )
