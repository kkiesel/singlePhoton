#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from sys import stdout

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
			stdout.write( "\r%s: %s / %s"%(fileName, event.GetReadEntry(), event.GetEntries() ) )
			stdout.flush()

		b = h_weight.FindBin( event.photons.at(0).ptJet(), event.ht )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()
	print

	f = ROOT.TFile( fileName, "update" )
	f.cd()
	weightTree.Write()
	f.Close()

def getMixedWeigthHisto( filenames, predFilenames, commonCut, control=True ):
	"""Calculate #photons/#photonFakes in bins of photons.ptJet and a second
	(global) variable.

	filenames: files containing photons
	predFilenames: files containing fakes
	"""

	regionCut = "met<100" if control else "met>=100"

	xVar = "photons[0].ptJet()"
	yVar = "ht"
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

	numerator = None
	for fileName in filenames:
		gTree = readTree( fileName, "photonTree" )
		num = createHistoFromTree2D( gTree, yVar+":"+xVar, "weight*( %s && %s )"%(regionCut, commonCut), xbinning, ybinning )
		if numerator:
			numerator.Add( num )
		else:
			numerator = num

	denominator = None
	for fileName in predFilenames:
		foTree = readTree( fileName, "photonJetTree" )
		den = createHistoFromTree2D( foTree, yVar+":"+xVar, "weight*( %s && %s )"%(regionCut, commonCut), xbinning, ybinning )
		if denominator:
			denominator.Add( den )
		else:
			denominator = den

	weight2D = divideHistos( numerator, denominator )

	# Set the weight and error for empty bins above the diagonal to one.
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if not weight2D.GetBinContent( i, j ) and weight2D.GetXaxis().GetBinLowEdge(i) < weight2D.GetYaxis().GetBinUpEdge(j):
				weight2D.SetBinContent( i, j, 1 )
				weight2D.SetBinError( i, j, 1 )

	drawWeightHisto( weight2D, numerator, denominator, control, shortName( filenames ) )
	return weight2D

def drawWeightHisto( weight2D, numerator, denominator, control, saveName ):
	regionString = "control" if control else "signal"
	# Draw the histograms
	info = PlotCaption(control=control, signal=not control)

	# Display the weight errors as 2D histograms.
	weight2D1 = divideHistos( numerator, denominator )
	weightErrors = weight2D.Clone( randomName() )
	weightErrors1 = weight2D.Clone( randomName() )
	weightRelErrors = weight2D.Clone( randomName() )
	weightRelErrors1 = weight2D.Clone( randomName() )
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if weight2D1.GetBinContent( i, j ):
				weightRelErrors.SetBinContent( i, j, weight2D.GetBinError( i, j )/weight2D.GetBinContent( i, j ) )
				weightRelErrors1.SetBinContent( i, j, weight2D1.GetBinError( i, j )/weight2D1.GetBinContent( i, j ) )

			weightErrors.SetBinContent( i, j, weight2D.GetBinError( i, j ) )
			weightErrors1.SetBinContent( i, j, weight2D1.GetBinError( i, j ) )

	# Draw histograms
	Styles.tdrStyle2D()
	ROOT.gStyle.SetPaintTextFormat("1.1f");
	can2D = ROOT.TCanvas()
	can2D.cd()
	for hist, name in [
			(numerator,"numerator"),
			(denominator, "denominator")]:
		hist.Draw("colz text")
		info.Draw()
		SaveAs(can2D, "qcd_preWeight_%s_%s_%s"%(saveName, name,regionString) )

	for hist, hist1, name in [
			(weight2D,weight2D1, "weight2D"),
			(weightErrors,weightErrors1, "weightError"),
			(weightRelErrors,weightRelErrors1, "weightRelError") ]:
		hist1.Draw("colz")
		hist.Draw("same text")
		info.Draw()
		SaveAs(can2D, "qcd_preWeight_%s_%s_%s"%(saveName,name,regionString) )

	Styles.tdrStyle()

	if control:
		weightFile = ROOT.TFile( "qcdWeight.root", "recreate" )
		weightFile.cd()
		weight2D.SetName("qcdWeight")
		weight2D.Write()
		weightErrors.Write()
		weightFile.Close()

def photonHisto( filenames, plot, cut, modifyEmptyBins ):
	gHist = None
	for filename in filenames:
		gTree = readTree( filename, "photonTree" )
		hist = getHisto( gTree, plot, cut=cut, color=1, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		if gHist:
			gHist.Add( hist )
		else:
			gHist = hist
	return gHist

def predictionHistos( filenames, plot, cut, modifyEmptyBins ):
	fHist, sysHist = None, None

	for filename in filenames:
		fTree = readTree( filename, "photonJetTree" )
		fTree.AddFriend( "foWeights", filename )

		# Prediction + statistical uncertainty coming from the loose control region
		hist = getHisto( fTree, plot, weight="weight*w_qcd", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )

		# This histogram contains the uncertainty due to the weight for small met
		# Here, the bin content is the uncertainty, the bin error is meaningless.
		sHist = getHisto( fTree, plot, weight="weight*w_qcd_error", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )

		# You can also do an approximation to calculate the error of the weight:
		doWeightApproximation = False
		if doWeightApproximation:
			sHistUp = getHisto( fTree, plot, weight="weight*(w_qcd_error+w_qcd)",
					cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
			sHistDown = getHisto( fTree, plot, weight="weight*(w_qcd-w_qcd_error)",
					cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
			for bin in range( sHistUp.GetNbinsX()+2 ):
				sHistUp.SetBinContent( bin, 0.5*(sHistUp.GetBinContent(bin)-sHistDown.GetBinContent(bin)))
			sHist = sHistUp
		# Approximation end ###################################################

		if fHist:
			fHist.Add( hist )
			sysHist.Add( sHist )
		else:
			fHist = hist
			sysHist = sHist

	sysHist.SetFillColor( sysHist.GetLineColor() )
	sysHist.SetLineColor( sysHist.GetLineColor() )
	sysHist.SetFillStyle(3254)
	sysHist.SetMarkerSize(0)

	return fHist, sysHist


def drawClosure( filenames, predFilenames, plot, commonCut, infoText, additionalLabel, modifyEmptyBins=True ):

	gHist = photonHisto( filenames, plot, commonCut, modifyEmptyBins )
	fHist, sysHist = predictionHistos( predFilenames, plot, commonCut, modifyEmptyBins )

	for bin in range( sysHist.GetNbinsX()+2 ):
		sysHist.SetBinError( bin, sysHist.GetBinContent(bin) )
		sysHist.SetBinContent( bin, fHist.GetBinContent(bin) )


	signalAbbrs = mergeDatasetAbbr( [ getDatasetAbbr(x) for x in filenames ] )
	predAbbrs =  mergeDatasetAbbr( [ getDatasetAbbr(x) for x in predFilenames ] )

	muhisto = Multihisto()
	muhisto.leg.SetHeader( ",".join( [ datasetToLatex(x) for x in signalAbbrs ] ) )
	muhisto.addHisto( gHist, "#gamma", draw="hist e0" )
	muhisto.addHisto( fHist, "#gamma_{jet}#upointw", draw="hist e0")
	muhisto.addHisto( sysHist, "#sigma_{w}", draw="e2")

	can = ROOT.TCanvas()
	can.cd()

	muhisto.Draw()

	fHistSumError = fHist.Clone( randomName() )
	for bin in range( fHistSumError.GetNbinsX()+2 ):
		fHistSumError.SetBinError( bin, sqrt( fHistSumError.GetBinError(bin)**2+sysHist.GetBinError(bin)**2 ) )
		gHist.SetBinError( bin, sqrt( gHist.GetBinError(bin)**2 + fHist.GetBinError(bin)**2 ) )

	from myRatio import Ratio
	#r = Ratio( "#gamma/#gamma_{pred}", gHist, fHistSumError )
	r = Ratio( "#gamma/#gamma_{pred}", gHist, sysHist )
	r.draw(0.75,1.25)
	infoText.Draw()
	saveLabelZeroBinError = "errorEmptyBin" if modifyEmptyBins else "normal"
	SaveAs(can, "qcdClosure_%s_%s"%("".join(signalAbbrs)+additionalLabel, plot) )

	# Since root is too stupid to clear the canvas before python is ending, clean
	# the canvas yourself
	ROOT.SetOwnership( can, False )
	del can



def qcdClosure( filenames, predFilenames, plots ):
	signalCut = "met>=100"
	controlCut = "!(%s)"%signalCut
	commonCut = "photons[0].ptJet()>110 && met < 100 && met > 50"
	#commonCut += " && photons[0].isGen(2)"
	#commonCut = "!photons[0].isGen(0) && !photons[0].isGen(1)"

	# Definition of labels
	infoControl = PlotCaption( control=True )
	infoSignal = PlotCaption( signal=True )
	info = PlotCaption()

	weights = getMixedWeigthHisto( filenames, predFilenames, commonCut )
	getMixedWeigthHisto( filenames, predFilenames, commonCut, False )

	for fileName in predFilenames:
		foTree = readTree( fileName, "photonJetTree" )
		writeWeight2DToFile( fileName, foTree, weights, "foWeights" )

	for plot in plots:
		if plot != "met":
			drawClosure( filenames, predFilenames, plot, commonCut+"&&"+signalCut, infoSignal, "_signal" )
			drawClosure( filenames, predFilenames, plot, commonCut+"&&"+controlCut, infoControl, "_control" )
		drawClosure( filenames, predFilenames, plot, commonCut, info, "", True )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("-p", "--prediction", nargs="+", default=[], type=isValidFile )
	arguments.add_argument("--plot", nargs="+", default=["met"] )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "photons[0].ptJet()","nGoodJets" ]
	if not opts.prediction:
		opts.prediction = opts.filenames

	qcdClosure( opts.filenames, opts.prediction, opts.plot )

