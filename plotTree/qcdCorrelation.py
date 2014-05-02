#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

def photonHisto( filenames, treename, plot, cut, modifyEmptyBins ):
	gHist = None
	for filename in filenames:
		gTree = readTree( filename, treename )
		hist = getHisto( gTree, plot, cut=cut, color=1, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		if gHist:
			gHist.Add( hist )
		else:
			gHist = hist
	return gHist


def qcdCorrelation( filenames, plot ):
	signalCut = "met>=100"
	controlCut = "!(%s)"%signalCut
	commonCut = " && photons[0].ptJet()>100 && photons[0].ptJet()<120 "

	gControlHist = photonHisto( filenames, "photonTree", plot, controlCut+commonCut, True )
	fControlHist = photonHisto( filenames, "photonJetTree", plot, controlCut+commonCut, True )
	gSignalHist = photonHisto( filenames, "photonTree", plot, signalCut+commonCut, True )
	fSignalHist = photonHisto( filenames, "photonJetTree", plot, signalCut+commonCut, True )
	gSignalHist.SetLineColor(2)
	fSignalHist.SetLineColor(2)
	fSignalHist.SetLineStyle(2)
	fControlHist.SetLineStyle(2)

	for h in [gControlHist, fControlHist, gSignalHist, fSignalHist]:
		h.Scale( 1./h.Integral() )
		h.SetMarkerSize(0)
		pass

	mh = Multihisto()
	mh.addHisto( gControlHist, "Control #gamma", draw="hist e" )
	mh.addHisto( fControlHist, "Control #gamma_{jet}",draw="hist e" )
	mh.addHisto( gSignalHist, "Signal #gamma",draw="hist e" )
	mh.addHisto( fSignalHist, "Signal #gamma_{jet}",draw="hist e" )



	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()
	info = PlotCaption()
	info.Draw()

	abbrs = mergeDatasetAbbr( [ getDatasetAbbr(x) for x in filenames ] )

	SaveAs( can, "correlation_%s_%s"%("".join(abbrs),plot ) )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--plot", nargs="+", default=["met"] )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "photons[0].ptJet()","Length$(jets.pt)", "Length$(photons.pt)"]

	for plot in opts.plot:
		qcdCorrelation( opts.filenames, plot )

