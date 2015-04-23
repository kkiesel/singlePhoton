#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multihisto import Multihisto
from treeFunctions import *
from myRatio import Ratio
from predictions import *

ROOT.gStyle.SetTitleSize( 0.035, "xyz" )
ROOT.gStyle.SetLabelSize( 0.035, "xyz" )



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
	#if plot != "met": commonCut += " && met>=100"

	gHist = getHists( filenames, plot, cut="(photons[0].isStatus(1)) && "+commonCut )
	eHist = multiDimFakeRate( filenames, plot, commonCut, isData=False )
	gHist.SetName("gHist")
	eHist.SetName("eHist")

	eHistSys = eHist.Clone( "eHistSys" )
	eHistSys = setRelativeUncertainty( eHistSys, 0.11 )
	eHistSys.SetFillColor( eHistSys.GetLineColor() )
	eHistSys.SetFillStyle(3354)
	eHistSys.SetMarkerSize(0)
	eHistSys.SetName( "eHistSys" )

	gDatasetAbbrs = [getDatasetAbbr(f) for f in filenames ]
	gDatasetAbbrs = mergeDatasetAbbr( gDatasetAbbrs )

	multihisto = Multihisto()
	multihisto.leg = myLegend(.63, .75, .96, .94 )
	multihisto.leg.SetTextSize(0.035)
	multihisto.leg.SetTextFont(42)

	multihisto.leg.SetHeader( ",".join([ datasetToLatex(x) for x in gDatasetAbbrs]) )
	multihisto.addHisto( gHist, "Direct Simulation", draw="p e x0" )
	multihisto.addHisto( eHist, "", draw="hist" )
	multihisto.addHisto( eHistSys, "Prediction", draw="e2" )

	infoText = ROOT.TLatex()
	infoText.SetNDC()
	infoText.SetTextFont( gHist.GetLabelFont() )
	infoText.SetTextSize( gHist.GetLabelSize() )
	infoText.SetText( .02, .96, "CMS Simulation                         19.7fb^{-1} (8 TeV)          #geq1#gamma,#geq2jets" )

	can = ROOT.TCanvas("c1", "", 600, 600 )
	can.cd()
	multihisto.Draw()
	infoText.Draw()
	r = Ratio( "Direct/Pred.", gHist, eHist, eHistSys )
	r.draw(0,2)
	SaveAs( can, "ewkClosure_%s_%s"%(getSaveNameFromDatasets(filenames), plot))
	ROOT.gPad.SaveAs( "plots/ewkClosure_%s_%s.C"%(getSaveNameFromDatasets(filenames), plot))

	ewkOut = ROOT.TFile("ewkOout.root", "recreate" )
	ewkOut.cd()
	for h in gHist, eHist, eHistSys:
		h.Write()
	ewkOut.Close()

	ROOT.SetOwnership( can, False )
	del can


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	opts = arguments.parse_args()

	closure( opts.filenames, opts.plot )

