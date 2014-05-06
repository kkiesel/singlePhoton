#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from myRatio import Ratio

def addRelativeUncertainty( hist, uncert ):
	for bin in range( hist.GetNbinsX()+1 ):
		hist.SetBinError( bin, hist.GetBinError(bin) |qPlus| (uncert * hist.GetBinContent(bin)) )
	return hist


def closure( filenames, predNames, errorZeroBins, opts ):
	commonCut = "!@muons.size() && !@electrons.size() && %s"%opts.cut
	for f in filenames:
		if "QCD_250" in f:
			commonCut += " && met>=100"

	gHist = None
	gDatasetAbbrs = []
	for filename in filenames:
		gDatasetAbbrs.append( getDatasetAbbr( filename ) )
		gTree = readTree( filename, "photonTree" )
		gamma = getHisto( gTree, opts.plot, color=1, cut="photons[0].isGen(1) && "+commonCut, fillEmptyBins=True )
		gamma.SetMarkerSize(0.3)

		if gHist:
			gHist.Add( gamma )
		else:
			gHist = gamma

	eHist = None
	eHistStat = None
	eDatasetAbbrs = []
	for filename in predNames:
		eDatasetAbbrs.append( getDatasetAbbr( filename ) )
		eTree = readTree( filename, "photonElectronTree" )

		recE = getHisto( eTree, opts.plot, color=2, fillEmptyBins=True, cut=commonCut )
		recE.SetFillColor( recE.GetLineColor() )
		recE.SetFillStyle(3354)
		recE.SetMarkerSize(0)
		recEStat = recE.Clone(randomName())
		recEStat = applyFakeRateEWK( recEStat, 0.0084, 0.000000001 )
		recEStat.SetFillStyle(0)
		recE = applyFakeRateEWK( recE, 0.0084, 0.0045 )

		if eHist:
			eHist.Add( recE )
			eHistStat.Add( recEStat )
		else:
			eHist = recE
			eHistStat = recEStat

	eHist = multiDimFakeRate( predNames )
	# needed if "setErrorInEmptyBins" enabled
	for h in [eHist, gHist, eHistStat]:
		bin = h.FindBin(99) # met bin before 100
		if not h.GetBinContent(bin):
			h.SetBinError( bin, 0 )

		h.GetXaxis().SetNdivisions(5,5,0, False)


	gDatasetAbbrs = mergeDatasetAbbr( gDatasetAbbrs )
	eDatasetAbbrs = mergeDatasetAbbr( eDatasetAbbrs )
	multihisto = Multihisto()
	if gDatasetAbbrs == eDatasetAbbrs:
		if len(gDatasetAbbrs) < 4:
			multihisto.leg.SetHeader( "/".join([ datasetToLatex(x) for x in gDatasetAbbrs ]) )
		else:
			multihisto.leg.SetHeader( "Total simulation" )
			multihisto.leg.SetX1( 0.6 )

		multihisto.addHisto( eHist, "Prediction", draw="e2" )
		multihisto.addHisto( gHist, "Simulation", draw="e0 hist" )
	else:
		multihisto.leg.SetX1(0.4)
		multihisto.addHisto( eHist, "#gamma_{e}#upointf^{,}_{e#rightarrow#gamma} "+",".join([ datasetToLatex(x) for x in eDatasetAbbrs ]), draw="e2" )
		multihisto.addHisto( gHist, "#gamma "+",".join([ datasetToLatex(x) for x in gDatasetAbbrs ]), draw="e0 hist" )



	infoText = ROOT.TLatex(0.03,.96, "CMS Private Work - 8TeV #geq1#gamma,#geq2jets" )
	infoText.SetNDC()

	can = ROOT.TCanvas()
	can.cd()
	multihisto.Draw()
	#eHistStat.Draw("same hist e")
	#gHist.Draw("same e0 hist")
	infoText.Draw()
	r = Ratio( "Sim./Pred.", gHist, eHist )
	r.draw(0,2)
	SaveAs( can, "ewkClosure_%s_%s_%s"%(getSaveNameFromDatasets(eDatasetAbbrs),opts.plot.replace(".",""),opts.save))



if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "-p", "--prediction", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--noError", action='store_false')
	arguments.add_argument( "--cut", default="1" )
	arguments.add_argument( "--save", default="new" )
	opts = arguments.parse_args()

	if not opts.prediction:
		opts.prediction = opts.filenames

	closure( opts.filenames, opts.prediction, opts.noError, opts )

