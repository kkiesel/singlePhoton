#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from myRatio import Ratio

def closure( filenames, predNames, errorZeroBins, opts ):
	gHist = None
	gDatasetAbbrs = []
	for filename in filenames:
		gDatasetAbbrs.append( getDatasetAbbr( filename ) )
		gTree = readTree( filename, "photonTree" )
		gamma = getHisto( gTree, opts.plot, color=1, weight="1",cut="photons.isGen(1) && @photons.size()>0", fillEmptyBins=True )

		if gHist:
			gHist.Add( gamma )
		else:
			gHist = gamma

	eHist = None
	eDatasetAbbrs = []
	for filename in predNames:
		eDatasetAbbrs.append( getDatasetAbbr( filename ) )
		eTree = readTree( filename, "photonElectronTree" )

		recE = getHisto( eTree, opts.plot, color=2, weight="1", fillEmptyBins=True )
		recE.SetFillColor( recE.GetLineColor() )
		recE.SetFillStyle(3254)
		recE.SetMarkerSize(0)
		recE = applyFakeRateEWK( recE )

		if eHist:
			eHist.Add( recE )
		else:
			eHist = recE

	gDatasetAbbrs = mergeDatasetAbbr( gDatasetAbbrs )
	eDatasetAbbrs = mergeDatasetAbbr( eDatasetAbbrs )
	multihisto = Multihisto()
	if gDatasetAbbrs == eDatasetAbbrs:
		multihisto.leg.SetHeader( "/".join([ datasetToLatex(x) for x in gDatasetAbbrs ]) )
		multihisto.addHisto( eHist, "#gamma_{e}#upointf'_{e#rightarrow#gamma}", draw="e2" )
		multihisto.addHisto( gHist, "#gamma", draw="e0 hist" )
	else:
		multihisto.leg.SetX1(0.4)
		multihisto.addHisto( eHist, "#gamma_{e}#upointf^{,}_{e#rightarrow#gamma} "+",".join([ datasetToLatex(x) for x in eDatasetAbbrs ]), draw="e2" )
		multihisto.addHisto( gHist, "#gamma "+",".join([ datasetToLatex(x) for x in gDatasetAbbrs ]), draw="e0 hist" )

	can = ROOT.TCanvas()
	hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
	hPad.cd()
	multihisto.Draw()

	ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
	ratioPad.cd()
	ratioPad.SetLogy(0)
	r = Ratio( "#gamma/#gamma_{pred}", gHist, eHist )
	ratio, sys, one = r.draw(0,2)
	ratio.Draw("same e0")
	sys.Draw("same e2")
	one.Draw()

	can.cd()
	hPad.Draw()
	ratioPad.Draw()
	SaveAs( can, "ewkClosure_%s_%s_%s"%("".join(eDatasetAbbrs),opts.plot.replace(".",""),opts.save))
	ROOT.SetOwnership(hPad, False)
	ROOT.SetOwnership(ratioPad, False)


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "-p", "--prediction", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--noError", action='store_false')
	arguments.add_argument( "--save", default="new" )
	opts = arguments.parse_args()

	if not opts.prediction:
		opts.prediction = opts.filenames

	closure( opts.filenames, opts.prediction, opts.noError, opts )

