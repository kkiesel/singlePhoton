#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from predictions import *

def getSignalHisto( scan="W", mg=1200, ms=1220, plot="met" ):
	if scan not in ["W", "B"]:
		print "Scan must be W or B"
		return
	if mg<0 or ms<0:
		print "Please insert sensible masses"

	plotRenaming = {"met": "gMet", "ht": "gHt", "nGoodJets": "gNjets", "photons[0].ptJet()": "gPt" }
	try:
		plotName = plotRenaming[plot]
	except LookupError:
		print "Signal Histogram not avaiable"
		return

	signalVersion="V03.00"
	h = readHisto( scan+"_gsq_%s.root"%signalVersion, "%s%s_%s"%(plotName,mg,ms) )

	# scale histo
	nGen = 60000 if scan == "W" else 10000
	w = getLumiWeight( "%s_%s_%s_375"%(scan,mg,ms), nGen, 19712 )
	h.Scale( w, "width" )

	# rebin
	import array
	binning = readAxisConf( plot )[2]
	xBins = array.array( 'd', binning )
	newh = h.Rebin( len(xBins)-1, randomName(), xBins )

	from inheritRoot import H1F
	newh.__class__ = H1F
	newh.MergeOverflow()

	return newh

def histogramToList( histo, error=False, minimalVal=100, roundInt=False ):
	# roundInt is necessary since ROOT has a low precission, leading to values
	# close to an integer (eg 23.000048)

	minimalBin = histo.FindBin(100)

	out = []

	# exclude under and overflow bin
	for bin in range(1, histo.GetNbinsX()+1):
		if bin < minimalBin:
			continue

		if error:
			out.append( histo.GetBinError(bin)*histo.GetBinWidth(bin) )
		else:
			if roundInt:
				out.append( int(round(histo.GetBinContent(bin)*histo.GetBinWidth(bin))) )
			else:
				out.append( histo.GetBinContent(bin)*histo.GetBinWidth(bin) )
	return out

def writeDataCard( versionData, dataHist, fgammaHist, fgammaWeightError,
		egammaHist, egammaHistsys,
		fsrZ, fsrZsys, fsrW, fsrWsys, fsrT, fsrTsys ):

	totalISR = addHistos( [fsrZ, fsrW, fsrT] )
	totalISRsys = addHistos( [fsrZsys, fsrWsys, fsrTsys] )

	# Print out events bin by bin
	binList = [100.0,]
	for bin in range( dataHist.GetNbinsX()+1 ):
		xval = dataHist.GetBinLowEdge(bin)
		if xval > 100:
			binList.append( xval )

	import time
	dataCardString = """
# This file contains event yield and uncertainties for data and simulated background event.
# Produced by:
#    Knut Kiesel
# on:
#    %s
# with input version:
#    %s
# Event yields so far are the sum of the background and not real data, since this
# analysis is still blind.
# For bin edges, the ROOT conversion is used (lower bin edge included, upper bin
# edge excluded).
"""%( time.strftime("%Y-%m-%d %H:%M:%S"), versionData )

	commonInformation = """
lumi = 19712 #pb
nMetBins = %s
""" % len(binList)

	binListInf = binList + ['infinity']

	for i in range(len(binList)):
		commonInformation += "bin %s = %s to %s\n"%(i,binListInf[i],binListInf[i+1])

	dataCardString += "\n##################################\n"
	dataCardString += commonInformation
	dataCardString += "\n##################################\n"

	data = histogramToList( dataHist, roundInt=True )
	qcd = histogramToList( fgammaHist )
	qcd_stat = histogramToList( fgammaHist, True )
	qcd_syst = histogramToList( fgammaWeightError, True )

	ewk = histogramToList( egammaHist )
	ewk_stat = histogramToList( egammaHist, True )
	ewk_syst = histogramToList( egammaHistsys, True )

	isr = histogramToList( totalISR )
	isr_stat = histogramToList( totalISR, True )
	isr_syst = histogramToList( totalISRsys, True )

	for name, array in [
		( "selected        = %s", data ),
		( "QCD background  = %s", qcd ),
		( "QCD stat uncert = %s", qcd_stat ),
		( "QCD syst uncert = %s", qcd_syst ),
		( "EWK background  = %s", ewk ),
		( "EWK stat uncert = %s", ewk_stat ),
		( "EWK syst uncert = %s", ewk_syst ),
		( "ISR background  = %s", isr ),
		( "ISR stat uncert = %s", isr_stat ),
		( "ISR syst uncert = %s", isr_syst ) ]:
		dataCardString += name%( ' '.join(map(str, array ) ) ) + "\n"

	#print dataCardString
	additionalInfo = "" if len(binList) == 6 else "_%smetBins"%len(binList)

	dataCardFileName = "eventYieldData%s-%s.txt"%(additionalInfo, time.strftime("%Y-%m-%d"))
	dataCardFile = open( dataCardFileName, "w")
	dataCardFile.write( dataCardString )
	dataCardFile.close()
	print "Write event yields to %s"%dataCardFileName

	return

def finalDistributionData( plot ):

	# Sample names
	mcVersion = "13"
	wg1 = "slimWGamma_50_130_V03.%s_tree.root"%mcVersion
	wg2 = "slimWGamma_130_inf_V03.%s_tree.root"%mcVersion
	tg = "slimTTGamma_V03.%s_tree.root"%mcVersion
	zgn = "slimZGammaNuNu_V03.%s_tree.root"%mcVersion

	published = False
	setRange = [ "A", "B" ] if published else [ "A", "B", "C", "D" ]
	additionalCut = "&& runNumber < 195948" if published else ""

	versionData = 13
	data = [ "PhotonHad%s_V03.%s_tree.root"%(x,versionData) for x in setRange ]

	leptonPtCut = 15 # only larger than 15 make sense here, since this is the reprocessing cut
	#commonCut = "(!@electrons.size() || Max$(electrons.pt)<{0}) && (!@muons.size() || Max$(muons.pt)<{0})".format(leptonPtCut)
	commonCut = "!@electrons.size() && !@muons.size()"
	#commonCut = "1"
	commonCut += additionalCut

	# Compute the weights:
	weight2D = getMixedWeigthHisto( data, data, commonCut )
	attachWeightsToFiles( data, weight2D, "foWeights" )
	from qcdClosure import drawWeightHisto
	drawWeightHisto( weight2D, "Data" )

	# Get Histograms
	dataHist = getHists( data, plot, commonCut+gCut )
	fgammaHist, fgammaWeightError = predictionHistos( data, plot, commonCut, modifyEmptyBins=False )

	egammaHist = multiDimFakeRate( data, plot, commonCut+gCut )

	fsrZ = getHists( [zgn], plot, commonCut+gCut )
	fsrW = getHists( [wg1,wg2], plot, commonCut+gCut )
	fsrT = getHists( [tg], plot, commonCut+gCut )

	#signal1 = getMetHisto( "W", 900, 1720 )
	#signal2 = getMetHisto( "B", 1700, 1120 )
	signal1 = getHists( ["slimW_1700_720_375_V03.06_tree.root"], plot, commonCut )
	signal2 = getHists( ["slimW_900_1720_375_V03.06_tree.root"], plot, commonCut )
	signal3 = getHists( ["slimB_1300_1720_375_V03.06_tree.root"], plot, commonCut )
	signal4 = getHists( ["slimB_1700_1120_375_V03.06_tree.root"], plot, commonCut )
	for i, signal in enumerate([signal1, signal2, signal3, signal4]):
		signal.SetLineColor( ROOT.kGreen +4 )
		signal.SetLineStyle(5+i)

	fgammaHist.SetLineColor(7)
	egammaHist.SetLineColor( 3 )
	fsrZ.SetLineColor( ROOT.kRed-7 )
	fsrW.SetLineColor( ROOT.kRed-9 )
	fsrT.SetLineColor( ROOT.kRed )

	mh = Multihisto()
	mh.orderByIntegral = False
	mh.addHisto( fsrT, "#gamma t#bar{t}", True )
	mh.addHisto( fsrW, "#gamma W", True )
	mh.addHisto( fsrZ, "#gamma Z", True )
	mh.addHisto( egammaHist, "e#rightarrow#gamma", True )
	mh.addHisto( fgammaHist, "Multijet", True )
	dataLegName = "Data"
	mh.addHisto( dataHist, dataLegName, draw="pe" )
	mh.addHisto( signal2, "Bino-like #chi_{1}^{0}", draw="hist" )
	mh.addHisto( signal1, "Wino-like #chi_{1}^{0}", draw="hist" )


	# additional ISR uncertainty
	ewkUncertainty = 0.11
	isrUncertaintyZ = 0.5
	isrUncertaintyW = 0.5
	isrUncertaintyT = 0.5

	# get all SYSTEMATICAL uncertainties:
	egammaHistsys = setRelativeUncertainty( egammaHist.Clone(randomName()), ewkUncertainty )
	fsrTsys = setRelativeUncertainty( fsrT.Clone(randomName()), isrUncertaintyT  )
	fsrWsys = setRelativeUncertainty( fsrW.Clone(randomName()), isrUncertaintyW  )
	fsrZsys = setRelativeUncertainty( fsrZ.Clone(randomName()), isrUncertaintyZ  )
	systematicUncertHistStack = ROOT.THStack()
	systematicUncertHistStack.Add( fgammaWeightError )
	systematicUncertHistStack.Add( egammaHistsys )
	systematicUncertHistStack.Add( fsrTsys )
	systematicUncertHistStack.Add( fsrWsys )
	systematicUncertHistStack.Add( fsrZsys )

	if plot == "met":
		writeDataCard( versionData, dataHist, fgammaHist, fgammaWeightError,
			egammaHist, egammaHistsys,
			fsrZ, fsrZsys, fsrW, fsrWsys, fsrT, fsrTsys )

	# draw stuff
	luminosity = 19.7
	infoText = ROOT.TLatex(0,.96, "CMS Private Work - %sfb^{-1} #sqrt{s}=8TeV #geq1#gamma_{tight},#geq2jets"%luminosity )
	infoText.SetNDC()
	infoText.SetTextSize(.04)

	can = ROOT.TCanvas()
	mh.Draw()
	errorBand = mh.stack.GetStack().Last().Clone( randomName() )
	sysUncert = systematicUncertHistStack.GetStack().Last()
	for bin in range(errorBand.GetNbinsX()+2):
		errorBand.SetBinError( bin, errorBand.GetBinError(bin) |qPlus| sysUncert.GetBinError(bin) )
	errorBand.SetFillStyle(3002)
	errorBand.SetMarkerSize(0)
	errorBand.SetFillColor(1)
	errorBand.Draw("e2 same")

	from myRatio import Ratio
	r = Ratio( "Data / Bkg", dataHist, mh.stack.GetStack().Last(), systematicUncertHistStack.GetStack().Last() )
	r.draw(0.5,1.5)

	infoText.Draw()

	SaveAs( can, "finalDistributionData_%s"%plot )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("--plot", nargs="+", default = ["met"] )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "photons[0].ptJet()","Length$(jets.pt)", "Length$(photons.pt)"]

	for plot in opts.plot:
		finalDistributionData( plot )

