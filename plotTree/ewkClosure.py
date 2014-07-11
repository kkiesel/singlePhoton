#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
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
	#commonCut ="1"

	gHist = getHists( filenames, plot, cut="(photons[0].isGen(1)||photons[0].isStatus(1)) && "+commonCut )
	eHist = multiDimFakeRate( filenames, plot, commonCut, isData=False )

	eHistSys = eHist.Clone( randomName() )
	eHistSys = setRelativeUncertainty( eHistSys, 0.11 )
	eHistSys.SetFillColor( eHistSys.GetLineColor() )
	eHistSys.SetFillStyle(3354)
	eHistSys.SetMarkerSize(0)

	gDatasetAbbrs = [getDatasetAbbr(f) for f in filenames ]
	gDatasetAbbrs = mergeDatasetAbbr( gDatasetAbbrs )

	multihisto = Multihisto()
	multihisto.leg = myLegend(.6, .75, .95, .92 )
	multihisto.leg.SetTextSize(0.035)
	multihisto.leg.SetTextFont(42)

	multihisto.leg.SetHeader( ",".join([ datasetToLatex(x) for x in gDatasetAbbrs]) )
	multihisto.addHisto( gHist, "Direct Simulation", draw="p e x0" )
	multihisto.addHisto( eHist, "Prediction", draw="hist" )
	multihisto.addHisto( eHistSys, "", draw="e2" )

	infoText = ROOT.TLatex()
	infoText.SetNDC()
	infoText.SetTextFont( gHist.GetLabelFont() )
	infoText.SetTextSize( gHist.GetLabelSize() )
	infoText.SetText( .02, .96, "CMS Simulation                       #sqrt{s}=8TeV, #intLdt=19.7fb^{-1}, #geq1#gamma,#geq2jets" )

	can = ROOT.TCanvas("", "", 600, 600 )
	can.cd()
	multihisto.Draw()
	infoText.Draw()
	r = Ratio( "Direct/Pred.", gHist, eHist, eHistSys )
	r.draw(0,2)
	SaveAs( can, "ewkClosure_%s_%s"%(getSaveNameFromDatasets(filenames), plot))

	gHist.SetName( "Direct Simulation" )
	eHist.SetName( "Prediction" )
	eHistSys.SetName( "weight" )
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

