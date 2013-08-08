#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
import argparse
import ConfigParser
import re
from math import sqrt
from multiplot import *
from treeFunctions import *
import ratios
import Styles

Styles.tdrStyle()
ROOT.gSystem.Load("libTreeObjects.so")

# to use user defined help message, sys.arv has to be sent to python and not
# to TApplication
ROOT.PyConfig.IgnoreCommandLineOptions = True

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
			print "%s / %s"%(event.GetReadEntry(), event.GetEntries() )
		b = h_weight.FindBin( event.photons.at(0).ptJet(), eval("event.%s"%reweightVar) )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()

	f = ROOT.TFile( fileName, "update" )
	f.cd()
	weightTree.Write()
	f.Close()

def getWeightHisto2D( gControlTree, foControlTree, datasetAffix ):
	"""The histogram for the weights in created here."""

	xVar = "photons[0].ptJet()"
	yVar = reweightVar
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

	weight_numerator = createHistoFromTree2D( gControlTree, yVar+":"+xVar, "weight*(met<100)", xbinning, ybinning )
	weight_denominator = createHistoFromTree2D( foControlTree, yVar+":"+xVar, "weight*(met<100)", xbinning, ybinning )
	weight2D = divideHistos( weight_numerator, weight_denominator )

	# Set the weight and error for empty bins to one.
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if not weight2D.GetBinContent( i, j ):
				weight2D.SetBinContent( i, j, 1 )
				weight2D.SetBinError( i, j, 1 )

	# Draw the histograms
	info = PlotCaption()
	info.controlCut()

	# Display the weight errors as 2D histograms.
	weightErrors = weight2D.Clone( randomName() )
	weightRelErrors = weight2D.Clone( randomName() )
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			weightErrors.SetBinContent( i, j, weight2D.GetBinError( i, j ) )
			if weight2D.GetBinContent( i, j ):
				weightRelErrors.SetBinContent( i, j, weight2D.GetBinError( i, j )/weight2D.GetBinContent( i, j ) )

	Styles.tdrStyle2D()
	ROOT.gStyle.SetPaintTextFormat("4.2f");
	can2D = ROOT.TCanvas()
	can2D.cd()
	for hist, name in [
			(weight_numerator,"numerator"),
			(weight_denominator, "denominator"),
			(weight2D, "weight2D"),
			(weightErrors, "weightError"),
			(weightRelErrors, "weightRelError") ]:
		hist.Draw("colz text")
		info.Draw()
		can2D.SaveAs( "plots/qcd_preWeight_%s_%s_%s.pdf"%(datasetAffix,name,yVar) )
	Styles.tdrStyle()

	return weight2D

def addErrorAndDivide( h_fo, h_fo_error ):
	import math
	h_new = h_fo_error.Clone( randomName() )
	for bin in range( h_fo_error.GetNbinsX()+1 ):
		if h_fo.GetBinContent(bin):
			h_new.SetBinError( bin, math.sqrt( h_fo_error.GetBinError(bin)**2 + h_fo.GetBinError(bin)**2 )/h_fo.GetBinContent(bin) )
		h_new.SetBinContent( bin, 1 )
	return h_new

def drawBeforeClosure( plot, gTree, foTree, cut, can, info, datasetAffix, additionalInfo="" ):

		# The first attempt to get the histogram is only to get the minimal
		# and maximal value on the x-axis, for not predefined binning
		h_gamma = getHisto( gTree, plot, cut=cut )
		h_fo = getHisto( foTree, plot, cut=cut )
		xMin, xMax = getXMinXMax( [ h_gamma, h_fo ] )

		h_gamma = getHisto( gTree, plot, cut=cut, color=1, firstBin=xMin,lastBin=xMax )
		h_fo = getHisto( foTree, plot, cut=cut, color=46, firstBin=xMin,lastBin=xMax )

		muhisto = Multihisto()
		muhisto.leg.SetHeader( datasetToLatex( datasetAffix ) )
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
		ratioGraph.hAxis.SetYTitle( "#gamma/#gamma_{pred}")

		can.cd()
		hPad.Draw()
		ratioPad.Draw()
		info.Draw()
		SaveAs(can, "plots","qcd_preWeight_%s_%s_%s"%(datasetAffix+additionalInfo, plot,reweightVar) )
		ROOT.SetOwnership( hPad, False )
		ROOT.SetOwnership( ratioPad, False )

def drawClosure( plot, gTree, foTree, cut, can, info, datasetAffix, additionalInfo="" ):

		# The first attempt to get the histogram is only to get the minimal
		# and maximal value on the x-axis, for not predefined binning
		h_gamma = getHisto( gTree, plot,cut=cut )
		h_fo = getHisto( foTree, plot, cut=cut )
		xMin, xMax = getXMinXMax( [ h_gamma, h_fo ] )

		h_gamma = getHisto( gTree, plot, cut=cut, color=1, firstBin=xMin,lastBin=xMax )
		h_fo = getHisto( foTree, plot, cut=cut, weight="weight*w_qcd", color=46, firstBin=xMin,lastBin=xMax )
		h_fo_error = getQCDErrorHisto( foTree, plot, cut=cut, firstBin=xMin, lastBin=xMax )
		h_fo_error.SetFillColor( h_fo.GetLineColor() )
		h_fo_error.SetLineColor( h_fo_error.GetLineColor() )
		h_fo_error.SetFillStyle(3254)
		h_fo_error.SetMarkerSize(0)

		muhisto = Multihisto()
		muhisto.leg.SetHeader( datasetToLatex( datasetAffix ) )
		muhisto.addHisto( h_gamma, "#gamma", draw="hist e0" )
		muhisto.addHisto( h_fo, "#gamma_{jet}#upointw", draw="hist e0")
		muhisto.addHisto( h_fo_error, "#sigma_{w}", draw="e2")

		hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
		hPad.cd()
		muhisto.Draw()

		ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
		ratioPad.cd()
		ratioPad.SetLogy( False )
		ratioGraph = ratios.RatioGraph( h_gamma, h_fo )
		ratioGraph.draw(ROOT.gPad, yMin=0, yMax=2, adaptiveBinning=False, errors="yx")
		ratioGraph.graph.Draw("same p e0") # draw nice points
		ratioGraph.hAxis.SetYTitle( "#gamma/#gamma_{pred}")

		h2 = addErrorAndDivide( h_fo, h_fo_error )
		h2.Draw("same e2")

		can.cd()
		hPad.Draw()
		ratioPad.Draw()
		info.Draw()
		SaveAs(can, "plots","qcd_afterWeighting_%s_%s_%s"%(datasetAffix+additionalInfo, plot,reweightVar) )
		ROOT.SetOwnership( hPad, False )
		ROOT.SetOwnership( ratioPad, False )

def qcdClosure( fileName, opts ):
	try:
		datasetAffix = re.match(".*slim(.*)_V.*", fileName ).groups()[-1]
	except:
		datasetAffix = fileName

	signalCut = "met>=100"
	controlCut = "!(%s)"%signalCut

	gTree = readTree( fileName, "photonTree" )
	foTree = readTree( fileName, "photonJetTree" )

	plots = [ "met" , "ht", "htHLT", "st80", "st30","photons[0].ptJet()","Length$(jets.pt)" ]

	# Definition of labels
	infoControl = PlotCaption()
	infoControl.controlCut()
	infoSignal = PlotCaption()
	infoSignal.signalCut()
	info = PlotCaption()

	can = ROOT.TCanvas()
	can.cd()

	for plot in plots:
		drawBeforeClosure( plot, gTree, foTree, signalCut, can, infoSignal, datasetAffix, "_signal" )
		drawBeforeClosure( plot, gTree, foTree, controlCut, can, infoControl, datasetAffix, "_control" )
		drawBeforeClosure( plot, gTree, foTree, "1", can, info, datasetAffix )

	commonCut = "htHLT > 500"
	gTree = gTree.CopyTree( commonCut )
	foTree = foTree.CopyTree( commonCut )

	weights = getWeightHisto2D( gTree, foTree, datasetAffix )
	writeWeight2DToFile( fileName, foTree, weights, "foWeights" )
	foTree.AddFriend( "foWeights", fileName )

	for plot in plots:
		drawClosure( plot, gTree, foTree, signalCut, can, infoSignal, datasetAffix, "_signal" )
		drawClosure( plot, gTree, foTree, controlCut, can, infoControl, datasetAffix, "_control" )
		drawClosure( plot, gTree, foTree, "1", can, info, datasetAffix )

	# Delete the trees from memory, important for many iterations.
	# The trees have to be deleted by hand to avoid segmentation violations.
	ROOT.SetOwnership( gTree, False )
	ROOT.SetOwnership( foTree, False )
	for tree in [ gTree, foTree ]:
		tree.Delete()
		del tree
	del can

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Calculate weighting "
			+"factors for QCD background estimation." )
	arguments.add_argument("--input", default=[""], nargs="+" )
	arguments.add_argument("--save", action="store_true", help="Save canvas as pdf.")
	opts = arguments.parse_args()

	for inName in opts.input:
		qcdClosure( inName, opts )

