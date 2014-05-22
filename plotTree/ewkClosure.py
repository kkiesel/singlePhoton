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

def closure( filenames, plot ):
	commonCut = "!@muons.size() && !@electrons.size()"
	commonCut ="1"

	gHist = getHists( filenames, plot, cut="photons[0].isGen(1) && "+commonCut )
	eHist = multiDimFakeRate( filenames, plot, commonCut, isData=False )

	eHistSys = eHist.Clone( randomName() )
	eHistSys = setRelativeUncertainty( eHistSys, 0.11 )
	eHistSys.SetFillColor( eHistSys.GetLineColor() )
	eHistSys.SetFillStyle(3354)
	eHistSys.SetMarkerSize(0)

	gDatasetAbbrs = [getDatasetAbbr(f) for f in filenames ]
	gDatasetAbbrs = mergeDatasetAbbr( gDatasetAbbrs )

	multihisto = Multihisto()
	multihisto.leg.SetHeader( "/".join([ datasetToLatex(x) for x in gDatasetAbbrs]) )
	multihisto.addHisto( gHist, "Simulation", draw="e0 hist" )
	multihisto.addHisto( eHist, "Prediction", draw="hist" )
	multihisto.addHisto( eHistSys, "", draw="e2" )

	infoText = ROOT.TLatex(0.03,.96, "CMS Private Work - 8TeV #geq1#gamma,#geq2jets" )
	infoText.SetNDC()

	can = ROOT.TCanvas()
	can.cd()
	multihisto.Draw()
	infoText.Draw()
	r = Ratio( "Sim./Pred.", gHist, eHist, eHistSys )
	r.draw(0,2)
	SaveAs( can, "ewkClosure_%s_%s"%(getSaveNameFromDatasets(filenames), plot))

	ROOT.SetOwnership( can, False )
	del can


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	opts = arguments.parse_args()

	closure( opts.filenames, opts.plot )

