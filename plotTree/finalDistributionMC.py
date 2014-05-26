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
		fsr, fsrsys ):

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

	isr = histogramToList( fsr )
	isr_stat = histogramToList( fsr, True )
	isr_syst = histogramToList( fsrsys, True )

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


def finalDistributionMC( plot ):

	# Sample names
	treeVersion = 12
	g1 = "slimGJets_200_400_V03.%s_tree.root"%treeVersion
	g2 = "slimGJets_400_inf_V03.%s_tree.root"%treeVersion
	q1 = "slimQCD_250_500_V03.%s_tree.root"%treeVersion
	q2 = "slimQCD_500_1000_V03.%s_tree.root"%treeVersion
	q3 = "slimQCD_1000_inf_V03.%s_tree.root"%treeVersion

	tt = "slimTTJets_V03.%s_tree.root"%treeVersion
	w1 = "slimWJets_250_300_V03.%s_tree.root"%treeVersion
	w2 = "slimWJets_300_400_V03.%s_tree.root"%treeVersion
	w3 = "slimWJets_400_inf_V03.%s_tree.root"%treeVersion


	wg1 = "slimWGamma_50_130_V03.%s_tree.root"%treeVersion
	wg2 = "slimWGamma_130_inf_V03.%s_tree.root"%treeVersion
	tg = "slimTTGamma_V03.%s_tree.root"%treeVersion
	zgn = "slimZGammaNuNu_V03.%s_tree.root"%treeVersion
	data = [ g1, g2, q1, q2, q3, tt, w1, w2, w3, zgn, wg1,wg2,tg ]

	# additional ISR uncertainty
	ewkUncertainty = 0.11
	isrUncertaintyZ = 0.5
	isrUncertaintyW = 0.5
	isrUncertaintyT = 0.5
	isrUncertainty = 0.5


	leptonPtCut = 25 # only larger than 15 make sense here, since this is the reprocessing cut
	commonCut = "(!@electrons.size() || Max$(electrons.pt)<{0}) && (!@muons.size() || Max$(muons.pt)<{0})".format(leptonPtCut)
	#commonCut = "!@electrons.size() && !@muons.size()"

	# Compute the weights:
	#weight2D = getMixedWeigthHisto( data, data, commonCut )
	#attachWeightsToFiles( data, weight2D, "foWeights" )
	#from qcdClosure import drawWeightHisto
	#drawWeightHisto( weight2D, "MC" )

	# Get Histograms
	dataHist = getHists( data, plot, commonCut )
	fgammaHist, fgammaWeightError = predictionHistos( data, plot, commonCut, modifyEmptyBins=False )

	egammaHist = multiDimFakeRate( data, plot, commonCut, isData=False )
	egammaHistsys = setRelativeUncertainty( egammaHist.Clone(randomName()), ewkUncertainty )

	fsrZ = getHists( [zgn], plot, commonCut )
	fsrZsys = setRelativeUncertainty( fsrZ, isrUncertaintyZ  )
	fsrW = getHists( [wg1,wg2], plot, commonCut )
	fsrWsys = setRelativeUncertainty( fsrW, isrUncertaintyW  )
	fsrT = getHists( [tg], plot, commonCut )
	fsrTsys = setRelativeUncertainty( fsrT, isrUncertaintyT  )
	fsr = addHistos( [fsrT, fsrW, fsrZ ] )
	fsrSys = setRelativeUncertainty( fsr, isrUncertainty )

	#signal1 = getMetHisto( "W", 900, 1720 )
	#signal2 = getMetHisto( "B", 1700, 1120 )
	signal1 = getHists( ["slimW_1700_720_375_V03.06_tree.root"], plot, commonCut )
	signal2 = getHists( ["slimW_900_1720_375_V03.06_tree.root"], plot, commonCut )
	signal3 = getHists( ["slimB_1300_1720_375_V03.06_tree.root"], plot, commonCut )
	signal4 = getHists( ["slimB_1700_1120_375_V03.06_tree.root"], plot, commonCut )
	for i, signal in enumerate([signal1, signal2, signal3, signal4]):
		signal.SetLineColor( ROOT.kGreen + i )
		signal.SetLineColor( ROOT.kBlue + i )
		#signal.SetLineStyle(5+i)

	# prettify histograms
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


	# get all SYSTEMATICAL uncertainties:
	systematicUncertHistStack = ROOT.THStack()
	systematicUncertHistStack.Add( fgammaWeightError )
	systematicUncertHistStack.Add( egammaHistsys )
	#systematicUncertHistStack.Add( fsrTsys )
	#systematicUncertHistStack.Add( fsrWsys )
	#systematicUncertHistStack.Add( fsrZsys )
	systematicUncertHistStack.Add( fsrSys )


	#if plot == "met":
	#	writeDataCard( treeVersion, dataHist, fgammaHist, fgammaWeightError,
	#		egammaHist, egammaHistsys,
	#		fsr, fsrSys )

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
		opts.plot = [ "met", "ht", "photons[0].ptJet()","nGoodJets", "@photons.size()"]

	for plot in opts.plot:
		finalDistributionMC( plot )

