#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from myRatio import Ratio
from predictions import *

def addRelativeUncertainty( hist, uncert ):
	for bin in range( hist.GetNbinsX()+1 ):
		hist.SetBinError( bin, hist.GetBinError(bin) |qPlus| (uncert * hist.GetBinContent(bin)) )
	return hist

def addRelativeUncertainty( hist, uncert ):
	for bin in range( hist.GetNbinsX()+1 ):
		hist.SetBinError( bin, uncert * hist.GetBinContent(bin) )
	return hist

def closure( filenames, plot, treename="photonJetTree" ):
	commonCut = "!@muons.size() && !@electrons.size()"
	#commonCut = "1"
	leptonPtCut = 15 # only larger than 15 make sense here, since this is the reprocessing cut
	commonCut = "(!@electrons.size() || Max$(electrons.pt)<{0}) && (!@muons.size() || Max$(muons.pt)<{0})".format(leptonPtCut)

	totalHist = getHists( filenames, plot, cut=commonCut, treeName=treename )
	gGenHist = getHists( filenames, plot, cut="photons[0].isGen(0) && "+commonCut, treeName=treename )
	eGenHist = getHists( filenames, plot, cut="photons[0].isGen(1) && "+commonCut, treeName=treename )
	gGenHist.SetLineColor(2)
	eGenHist.SetLineColor(3)

	gDatasetAbbrs = [getDatasetAbbr(f) for f in filenames ]
	gDatasetAbbrs = mergeDatasetAbbr( gDatasetAbbrs )

	multihisto = Multihisto()
	multihisto.leg.SetHeader( "/".join([ datasetToLatex(x) for x in gDatasetAbbrs]) )
	multihisto.addHisto( totalHist, "Simulation", draw="e0 hist" )
	multihisto.addHisto( gGenHist, "gen #gamma", toStack=True, draw="e0 hist" )
	multihisto.addHisto( eGenHist, "gen e", toStack=True, draw="e0 hist" )

	infoText = ROOT.TLatex(0.03,.96, "CMS Private Work - 8TeV #geq1#gamma,#geq2jets" )
	infoText.SetNDC()

	can = ROOT.TCanvas()
	can.cd()
	multihisto.Draw()
	infoText.Draw()
	r = Ratio( "Sim./Pred.", totalHist, multihisto.stack.GetStack().Last() )
	r.draw(0,2)
	SaveAs( can, "genMatches_%s_%s"%(getSaveNameFromDatasets(filenames), plot))

	ROOT.SetOwnership( can, False )
	del can


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	opts = arguments.parse_args()

	closure( opts.filenames, opts.plot )

