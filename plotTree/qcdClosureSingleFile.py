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
			stdout.write( "\r%s: %s / %s"%(fileName, event.GetReadEntry(), event.GetEntries() ) )
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

def getMixedWeigthHisto( filenames, predFilenames, commonCut ):
	"""Calculate #photons/#photonFakes in bins of photons.ptJet and a second
	(global) variable.

	filenames: files containing photons
	predFilenames: files containing fakes
	"""

	xVar = "photons[0].ptJet()"
	yVar = reweightVar
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

	numerator = None
	for fileName in filenames:
		gTree = readTree( fileName, "photonTree" )
		num = createHistoFromTree2D( gTree, yVar+":"+xVar, "weight*(met<100&&%s)"%commonCut, xbinning, ybinning )
		if numerator:
			numerator.Add( num )
		else:
			numerator = num

	denominator = None
	for fileName in predFilenames:
		foTree = readTree( fileName, "photonJetTree" )
		den = createHistoFromTree2D( foTree, yVar+":"+xVar, "weight*(met<100&&%s)"%commonCut, xbinning, ybinning )
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

	drawWeightHisto( weight2D, numerator, denominator )
	return weight2D

def drawWeightHisto( weight2D, numerator, denominator ):
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
	ROOT.gStyle.SetPaintTextFormat("1.1f");
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
		SaveAs(can2D, "qcd_preWeight_%s_%s"%(name,reweightVar) )
	Styles.tdrStyle()

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
		hist = getHisto( fTree, plot, weight="weight*w_qcd", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		histUp = getHisto( fTree, plot, weight="weight*(w_qcd+w_qcd_error)", cut=cut, color=hist.GetLineColor(), firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		histDown = getHisto( fTree, plot, weight="weight*(w_qcd-w_qcd_error)", cut=cut, color=hist.GetLineColor(), firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		sHist = histUp.Clone( randomName() )
		for bin in range( histUp.GetNbinsX()+2 ):
			up = histUp.GetBinContent(bin)
			down = histDown.GetBinContent(bin)
			sHist.SetBinContent( bin, (up+down)/2 )
			sHist.SetBinError( bin, (up-down)/2 )

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


def drawClosure( filenames, predFilenames, plot, commonCut, can, infoText, additionalLabel, modifyEmptyBins=False ):

	gHist = photonHisto( filenames, plot, commonCut, modifyEmptyBins )
	fHist, sysHist = predictionHistos( predFilenames, plot, commonCut, modifyEmptyBins )
	fewkHist, sysewkHist = predictionHistos( predFilenames+["slimTTJets_V02.22_tree.root", "slimWJets_V02.22_tree.root"], plot, commonCut, modifyEmptyBins )
	for h in [ fewkHist, sysewkHist]:
		h.SetLineColor(3)
		h.SetMarkerColor(3)
	sysewkHist.SetFillColor(3)
	sysewkHist.SetFillStyle(3254)
	sysewkHist.SetMarkerSize(0)

	fDatasetAbbrs = []
	for f in predFilenames:
		fDatasetAbbrs.append( getDatasetAbbr( f ) )


	fullDatasetAbbr = mergeDatasetAbbr( fDatasetAbbrs )

	muhisto = Multihisto()
	muhisto.leg.SetHeader( ",".join( [ datasetToLatex(x) for x in fullDatasetAbbr ] ) )
	muhisto.addHisto( gHist, "#gamma", draw="hist e0" )
	muhisto.addHisto( fHist, "#gamma_{jet}#upointw", draw="hist e0")
	muhisto.addHisto( sysHist, "#sigma_{w}", draw="e2")
	muhisto.addHisto( fewkHist, "pred with ewk", draw="hist e0")
	muhisto.addHisto( fewkHist, "pred with ewk", draw="e2")

	hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
	hPad.cd()
	muhisto.Draw()

	ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
	ratioPad.cd()
	ratioPad.SetLogy( False )
	from myRatio import Ratio
	r = Ratio( "#gamma/#gamma_{pred ewk}", gHist, fewkHist )
	ratio, sys, one = r.draw(0,2)
	ratio.Draw("same e0")
	sys.Draw("same e2")
	one.Draw()

	can.cd()
	hPad.Draw()
	ratioPad.Draw()
	infoText.Draw()
	saveLabelZeroBinError = "errorEmptyBin" if modifyEmptyBins else "normal"
	SaveAs(can, "qcdClosure_%s_%s_%s_%s"%("".join( fullDatasetAbbr)+additionalLabel, plot, saveLabelZeroBinError, reweightVar) )
	ROOT.SetOwnership( hPad, False )
	ROOT.SetOwnership( ratioPad, False )



def qcdClosure( filenames, predFilenames, plots ):
	signalCut = "met>=100"
	controlCut = "!(%s)"%signalCut
	commonCut = "!photons.isGen(1) && @photons.size()"
	commonCut = "1"

	# Definition of labels
	infoControl = PlotCaption( control=True )
	infoSignal = PlotCaption( signal=True )
	info = PlotCaption()

	can = ROOT.TCanvas()
	can.cd()

	weights = getMixedWeigthHisto( filenames, predFilenames, commonCut )

	for fileName in predFilenames:
		foTree = readTree( fileName, "photonJetTree" )
		writeWeight2DToFile( fileName, foTree, weights, "foWeights" )

	for plot in plots:
		if plot != "met":
			drawClosure( filenames, predFilenames, plot, commonCut+"&&"+signalCut, can, infoSignal, "_signal" )
			drawClosure( filenames, predFilenames, plot, commonCut+"&&"+controlCut, can, infoControl, "_control" )
		drawClosure( filenames, predFilenames, plot, commonCut, can, info,  "" )
		drawClosure( filenames, predFilenames, plot, commonCut, can, info, "", True )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("-p", "--prediction", nargs="+", default=[], type=isValidFile )
	arguments.add_argument("--plot", nargs="+", default=["met"] )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "photons[0].ptJet()","Length$(jets.pt)", "Length$(photons.pt)"]
	if not opts.prediction:
		opts.prediction = opts.filenames

	qcdClosure( opts.filenames, opts.prediction, opts.plot )

