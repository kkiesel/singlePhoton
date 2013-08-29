#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from sys import stdout

# variable to reweight
reweightVar = "ht"

def writeWeight2DToFile( fileName, tree, h_weight, weightTreeName ):
	"""Write weight for a tree in another tree in a given file.
	This tree can be added to the original tree via 'AddFriend()'.
	fileName: name of file to which tree is written
	tree: tree which is weighted
	h_weight: two dimensional histogram with weights
	weighTreeName: name of the new tree
	"""
	weightTree = ROOT.TTree( weightTreeName, "Tree containing QCD weights" )
	import numpy
	weight = numpy.zeros( 1, dtype=float)
	weight_error = numpy.zeros( 1, dtype=float)
	weightTree.Branch( "w_qcd", weight, "w_qcd/D" )
	weightTree.Branch( "w_qcd_error", weight_error, "w_qcd_error/D" )

	for event in tree:
		if not event.GetReadEntry()%10000:
			stdout.write( "\r%s / %s"%(event.GetReadEntry(), event.GetEntries() ) )
			stdout.flush()

		b = h_weight.FindBin( event.photons.at(0).ptJet(), eval("event.%s"%reweightVar) )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()
	print

	f = ROOT.TFile( fileName, "update" )
	f.cd()
	weightTree.Write()
	f.Close()

def drawBeforeClosure( plot, gTree, foTree, cut, can, info, datasetAbbr, additionalInfo="" ):
		# The first attempt to get the histogram is only to get the minimal
		# and maximal value on the x-axis, for not predefined binning
		h_gamma = getHisto( gTree, plot, cut=cut )
		h_fo = getHisto( foTree, plot, cut=cut )
		xMin, xMax = getXMinXMax( [ h_gamma, h_fo ] )

		h_gamma = getHisto( gTree, plot, cut=cut, color=1, firstBin=xMin,lastBin=xMax )
		h_fo = getHisto( foTree, plot, cut=cut, color=46, firstBin=xMin,lastBin=xMax )

		muhisto = Multihisto()
		muhisto.leg.SetHeader( datasetToLatex( datasetAbbr ) )
		muhisto.addHisto( h_gamma, "#gamma", draw="hist e0" )
		muhisto.addHisto( h_fo, "#gamma_{jet}", draw="hist e0")

		hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
		hPad.cd()
		muhisto.Draw()

		ROOT.TGaxis.SetMaxDigits(4)
		ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
		ratioPad.cd()
		ratioPad.SetLogy( False )
		ratioGraph = ratios.RatioGraph( h_gamma, h_fo )
		ratioGraph.draw(ROOT.gPad, yMin=None, yMax=None, adaptiveBinning=False, errors="yx")
		ratioGraph.graph.Draw("same p e0") # draw nice points
		ratioGraph.hAxis.SetYTitle( "#gamma/#gamma_{jet}")

		can.cd()
		hPad.Draw()
		ratioPad.Draw()
		info.Draw()
		SaveAs(can, "qcd_preWeight_%s_%s_%s"%(datasetAbbr+additionalInfo, plot,reweightVar) )
		ROOT.SetOwnership( hPad, False )
		ROOT.SetOwnership( ratioPad, False )

def getClosureHists( plot, gTree, foTree, cut, can, infoSignal ):
	gHist = getHisto( gTree, plot, cut=cut, color=1, firstBin=1, lastBin=1 )
	foHist = getHisto( foTree, plot, weight="weight*w_qcd", cut=cut, color=46, firstBin=1, lastBin=1 )
	upHist = getHisto( foTree, plot, weight="weight*(w_qcd+w_qcd_error)", cut=cut, color=foHist.GetLineColor(), firstBin=1, lastBin=1 )
	downHist = getHisto( foTree, plot, weight="weight*(w_qcd-w_qcd_error)", cut=cut, color=foHist.GetLineColor(), firstBin=1, lastBin=1 )

	return gHist, foHist, upHist, downHist

def getMixedWeigthHisto( filenames, commonCut, datasetAbbr ):

	xVar = "photons[0].ptJet()"
	yVar = reweightVar
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

	numerator = None
	denominator = None

	for fileName in filenames:
		gTree = readTree( fileName, "photonTree" )
		foTree = readTree( fileName, "photonJetTree" )

		num = createHistoFromTree2D( gTree, yVar+":"+xVar, "weight*(met<100&&%s)"%commonCut, xbinning, ybinning )
		den = createHistoFromTree2D( foTree, yVar+":"+xVar, "weight*(met<100&&%s)"%commonCut, xbinning, ybinning )

		if numerator:
			numerator.Add( num )
		else:
			numerator = num
		if denominator:
			denominator.Add( den )
		else:
			denominator = den

	weight2D = divideHistos( numerator, denominator )

	# Set the weight and error for empty bins to one.
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if not weight2D.GetBinContent( i, j ):
				weight2D.SetBinContent( i, j, 1 )
				weight2D.SetBinError( i, j, 1 )

	# Draw the histograms
	info = PlotCaption(control=True)

	# Display the weight errors as 2D histograms.
	weightErrors = weight2D.Clone( randomName() )
	weightRelErrors = weight2D.Clone( randomName() )
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			weightErrors.SetBinContent( i, j, weight2D.GetBinError( i, j ) )
			if weight2D.GetBinContent( i, j ):
				weightRelErrors.SetBinContent( i, j, weight2D.GetBinError( i, j )/weight2D.GetBinContent( i, j ) )

	# Draw histograms
	Styles.tdrStyle2D()
	ROOT.gStyle.SetPaintTextFormat("4.2f");
	can2D = ROOT.TCanvas()
	can2D.cd()
	for hist, name in [
			(numerator,"numerator"),
			(denominator, "denominator"),
			(weight2D, "weight2D"),
			(weightErrors, "weightError"),
			(weightRelErrors, "weightRelError") ]:
		hist.Draw("colz text")
		info.Draw()
		SaveAs(can2D, "qcd_preWeight_%s_%s_%s"%(datasetAbbr,name,yVar) )
	Styles.tdrStyle()

	weightFile = ROOT.TFile( "qcdWeight.root", "recreate" )
	weightFile.cd()
	weight2D.SetName("qcdWeight")
	weight2D.Write()
	weightErrors.Write()
	weightFile.Close()


	return weight2D

def drawClosure( filenames, plot, commonCut, can, infoText, fullDatasetAbbr, additionalLabel, errorZeroBins=False ):
	hists = {}
	for fileName in filenames:
		gTree = readTree( fileName, "photonTree" )
		foTree = readTree( fileName, "photonJetTree" )
		foTree.AddFriend( "foWeights", fileName )

		allHists = getClosureHists( plot, gTree, foTree, commonCut, can, infoText )

		#assume the mean weight is the same in gTree and foTree
		weightH = createHistoFromTree( gTree, "weight", "weight", 100 )
		weight = weightH.GetMean()

		poissonZeroError = 1.14 if errorZeroBins else 0

		for h in allHists:
			for bin in range(1, h.GetNbinsX()+2):
				# if the bin left or right is not empty but the bin itself, set the error
				if not h.GetBinContent( bin ) and ( h.GetBinContent( bin-1 ) or h.GetBinContent( bin+1 ) ):
					h.SetBinError( bin, poissonZeroError*weight / h.GetBinWidth(bin) )
		hists[ fileName ] = allHists

	addedList = [None, None, None, None]
	for name, hList in hists.iteritems():
		for i, h in enumerate(addedList):
			if h:
				addedList[i].Add( hList[i] )
			else:
				addedList[i] = hList[i]
	gHist, foHist, upHist, downHist = addedList

	sysHist = upHist.Clone( randomName() )
	for bin in range( upHist.GetNbinsX()+2 ):
		up = upHist.GetBinContent(bin)
		down = downHist.GetBinContent(bin)
		sysHist.SetBinContent( bin, (up+down)/2 )
		sysHist.SetBinError( bin, (up-down)/2 )

	sysHist.SetFillColor( foHist.GetLineColor() )
	sysHist.SetLineColor( sysHist.GetLineColor() )
	sysHist.SetFillStyle(3254)
	sysHist.SetMarkerSize(0)

	muhisto = Multihisto()
	muhisto.leg.SetHeader( datasetToLatex( fullDatasetAbbr ) )
	muhisto.addHisto( gHist, "#gamma", draw="hist e0" )
	muhisto.addHisto( foHist, "#gamma_{jet}#upointw", draw="hist e0")
	muhisto.addHisto( sysHist, "#sigma_{w}", draw="e2")

	hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
	hPad.cd()
	muhisto.Draw()

	ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
	ratioPad.cd()
	ratioPad.SetLogy( False )
	from myRatio import Ratio
	r = Ratio( "#gamma/#gamma_{pred}", gHist, foHist )
	ratio, sys, one = r.draw(0,2)
	ratio.Draw("same e0")
	sys.Draw("same e2")
	one.Draw()

	can.cd()
	hPad.Draw()
	ratioPad.Draw()
	infoText.Draw()
	SaveAs(can, "qcdClosure_%s_%s_%s_%s"%(fullDatasetAbbr+additionalLabel, plot,errorZeroBins, reweightVar) )


def qcdClosure( filenames, plots, fullDatasetAbbr ):
	signalCut = "met>=100"
	controlCut = "!(%s)"%signalCut
	commonCut = "1"

	# Definition of labels
	infoControl = PlotCaption( control=True )
	infoSignal = PlotCaption( signal=True )
	info = PlotCaption()

	can = ROOT.TCanvas()
	can.cd()

	weights = getMixedWeigthHisto( filenames, commonCut, fullDatasetAbbr )

	for fileName in filenames:
		foTree = readTree( fileName, "photonJetTree" )
		writeWeight2DToFile( fileName, foTree, weights, "foWeights" )

	for plot in plots:
		if plot != "met":
			drawClosure( filenames, plot, commonCut+"&&"+signalCut, can, infoSignal, fullDatasetAbbr, "_signal" )
			drawClosure( filenames, plot, commonCut+"&&"+controlCut, can, infoControl, fullDatasetAbbr, "_control" )
		drawClosure( filenames, plot, commonCut, can, info, fullDatasetAbbr, "" )
		drawClosure( filenames, plot, commonCut, can, info, fullDatasetAbbr, "", True )

	#ROOT.SetOwnership( hPad, False )
	#ROOT.SetOwnership( ratioPad, False )

	# Delete the trees from memory, important for many iterations.
	# The trees have to be deleted by hand to avoid segmentation violations.
	#ROOT.SetOwnership( gTree, False )
	#ROOT.SetOwnership( foTree, False )
	#for tree in [ gTree, foTree ]:
	#	tree.Delete()
	#	del tree
	del can

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--plot", nargs="+", default = ["met"] )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "photons[0].ptJet()","Length$(jets.pt)", "Length$(photons.pt)"]

	datasets = []
	for filename in opts.filenames:
		datasets.append( getDatasetAbbr( filename ) )
	datasets.sort()

	if "GJets_200_400" in datasets and "GJets_400_inf" in datasets:
		datasets.remove( "GJets_200_400" )
		datasets.remove( "GJets_400_inf" )
		datasets.append( "GJets" )
	if "QCD_250_500" in datasets and "QCD_500_1000" in datasets and "QCD_1000_inf" in datasets:
		datasets.remove( "QCD_250_500" )
		datasets.remove( "QCD_500_1000" )
		datasets.remove( "QCD_1000_inf" )
		datasets.append( "QCD" )
	if "QCD" in datasets and "GJets" in datasets:
		datasets.remove( "QCD" )
		datasets.remove( "GJets" )
		datasets.append( "AllQCD" )

	sumDatasetAbbr = "-".join( datasets )

	qcdClosure( opts.filenames, opts.plot, sumDatasetAbbr )

