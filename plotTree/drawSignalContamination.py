#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

def applyFakeRate( histo, f, e_f ):
	for i in range( histo.GetNbinsX() +1 ):
		binContent = histo.GetBinContent(i)
		histo.SetBinContent( i, histo.GetBinContent(i)*f)
		# sigma_{ef} = sqrt( ( e*e_f )**2 + ( e_e*f )**2 )
		histo.SetBinError( i, sqrt((binContent*e_f)**2 + (histo.GetBinError(i)*f)**2 ) )
	return histo

def drawSignalContamination( filename, plot ):

	plot = "met"

	###########################################################
	fakeRate = 0.0084
	fakeRateError = 0.0006

	# correct fake rate, if it is estimated with yutaros method
	fakeRateError = fakeRateError / (1-fakeRate)**2
	fakeRate = fakeRate / ( 1 - fakeRate )

	weight2D = readHisto( "qcdWeight.root", "qcdWeight" )
	###########################################################

	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(False)

	gTree = readTree( filename, "photonTree" )
	gHisto = getHisto( gTree, plot )

	eTree = readTree( filename, "photonElectronTree" )
	eHisto = applyFakeRate( getHisto( eTree, plot, color=2 ), fakeRate, fakeRateError )
	eContamination = divideHistos( eHisto, gHisto )
	eContamination.GetYaxis().SetTitle("signal contamination")

	fTree = readTree( filename, "photonJetTree" )
	from qcdClosureSingleFile import writeWeight2DToFile
	writeWeight2DToFile( filename, fTree, weight2D, "foWeights" )
	fTree.AddFriend( "foWeights", filename )
	fHisto = getHisto( fTree, plot, color=4, weight="weight*w_qcd" )
	fContamination = divideHistos( fHisto, gHisto )

	totalContamination = divideHistos( addHistos( [eHisto, fHisto] ), gHisto )
	totalContamination.SetLineColor(1)
	totalContamination.SetMarkerColor(1)

	mh = Multihisto()
	mh.leg.SetHeader("Prediction from")
	mh.addHisto( eContamination, "#gamma_{e}", draw="e0" )
	mh.addHisto( fContamination, "#gamma_{jet}", draw="e0" )
	mh.addHisto( totalContamination, "total", draw="e0" )
	mh.Draw()

	SaveAs(can, "signalContamination_%s"%getDatasetAbbr( filename ) )



if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--plot", nargs="+", default = ["met"] )
	opts = arguments.parse_args()

	for inName in opts.filenames:
		drawSignalContamination( inName, opts.plot )

