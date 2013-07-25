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
		b = h_weight.FindBin( event.photon.at(0).ptJet, event.ht )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()

	f = ROOT.TFile( fileName, "update" )
	f.cd()
	weightTree.Write()
	f.Close()

def drawDifferentWeights1D( weight2D ):
	"""Draw the weights as one dimensional histograms for each y-Bin."""
	weightHistos = []
	for yBin in range( weight2D.GetNbinsY()+1 ):
		h = weight2D.ProjectionX( "_px%s"%yBin, yBin, yBin+1 )
		h.SetLineColor( yBin+1 )
		h.SetMarkerColor( h.GetLineColor() )
		h.SetTitle("w")
		h.SetMarkerSize(0)
		weightHistos.append( h )

	can = ROOT.TCanvas()
	can.SetLogy(0)
	can.cd()

	weightHistos[0].Draw("hist e" )
	weightHistos[0].SetMaximum(33)
	for h in weightHistos[1:]:
		h.Draw("same hist e")

	can.SaveAs("qcd_preWeight_weights1D.pdf")
	return

def getWeightHisto2D( gControlTree, foControlTree, datasetAffix ):
	"""The histogram for the weights in created here."""

	xVar = "photon[0].ptJet"
	yVar = "ht"
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

	weight_numerator = createHistoFromTree2D( gControlTree, yVar+":"+xVar, "weight", xbinning, ybinning )
	weight_denominator = createHistoFromTree2D( foControlTree, yVar+":"+xVar, "weight", xbinning, ybinning )
	weight2D = divideHistos( weight_numerator, weight_denominator )

	# Set the weight and error for empty bins to one.
	"""for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if not weight2D.GetBinContent( i, j ):
				weight2D.SetBinContent( i, j, 1 )
				weight2D.SetBinError( i, j, 1 )
	"""

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
		can2D.SaveAs( "plots/qcd_preWeight_%s_%s.pdf"%(datasetAffix,name) )
	Styles.tdrStyle()

	# Draw the projections before the weights are applied.
	pt = getHisto( gControlTree, xVar )
	pt_fo = getHisto( foControlTree, xVar, color=2 )
	ht = getHisto( gControlTree, yVar )
	ht_fo = getHisto( foControlTree, yVar, color=2 )

	can = ROOT.TCanvas()
	can.cd()
	x_multi = Multihisto()
	x_multi.addHisto( pt, "#gamma", draw="hist" )
	x_multi.addHisto( pt_fo, "#gamma_{jet}", draw="hist" )
	x_multi.Draw()
	info.Draw()
	can.SaveAs( "plots/qcd_preWeight_%s_pt.pdf"%datasetAffix )

	y_multi = Multihisto()
	y_multi.addHisto( ht, "#gamma", draw="hist" )
	y_multi.addHisto( ht_fo, "#gamma_{jet}", draw="hist" )
	y_multi.Draw()
	info.Draw()
	can.SaveAs( "plots/qcd_preWeight_%s_ht.pdf"%datasetAffix )

	return weight2D

def clearTreeFromJets( inTree, jets, _deltaR=.3 ):
	from splitCandidatesQCDtest import deltaR
	newTree = inTree.CloneTree(0)
	newTree.SetBranchAddress( "jet", jets )
	for event in inTree:
		jets.clear()
		for j in event.jet:
			if deltaR( j, event.photon[0] ) > _deltaR:
				jets.push_back( j )
		if jets.size() > 1:
			newTree.Fill()

	# Set the size of the tree in memory to 10MB or get memory overflow.
	#newTree.SetMaxVirtualSize(int(1e7))
	return newTree

def addErrorAndDivide( h_fo, h_fo_error ):
	import math
	h_new = h_fo_error.Clone( randomName() )
	for bin in range( h_fo_error.GetNbinsX()+1 ):
		if h_fo.GetBinContent(bin):
			h_new.SetBinError( bin, math.sqrt( h_fo_error.GetBinError(bin)**2 + h_fo.GetBinError(bin)**2 )/h_fo.GetBinContent(bin) )
		h_new.SetBinContent( bin, 1 )
	return h_new

def drawBeforeClosure( plot, gSignalTree, foSignalTree, can, info, datasetAffix, additionalInfo="" ):

		# The first attempt to get the histogram is only to get the minimal
		# and maximal value on the x-axis, for not predefined binning
		h_gamma = getHisto( gSignalTree, plot )
		h_fo = getHisto( foSignalTree, plot, cut="1" )
		xMin, xMax = getXMinXMax( [ h_gamma, h_fo ] )

		h_gamma = getHisto( gSignalTree, plot, color=1, firstBin=xMin,lastBin=xMax )
		h_fo = getHisto( foSignalTree, plot, cut="1", color=46, firstBin=xMin,lastBin=xMax )

		if plot == "met" and additionalInfo=="":
			print "fo, met < 100: ", h_fo.Integral(0,h_fo.FindBin(100),"width")
			print "gamma, met < 100: ", h_gamma.Integral(0,h_gamma.FindBin(100),"width")
			print "fo, met > 100: ", h_fo.Integral(h_fo.FindBin(100),h_fo.GetNbinsX()+1,"width")
			print "gamma, met > 100: ", h_gamma.Integral(h_gamma.FindBin(100),h_gamma.GetNbinsX()+1,"width")


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
		SaveAs(can, "plots","qcd_preWeight_%s_%s"%(datasetAffix+additionalInfo, plot) )
		ROOT.SetOwnership( hPad, False )
		ROOT.SetOwnership( ratioPad, False )

def drawClosure( plot, gSignalTree, foSignalTree, can, info, datasetAffix, additionalInfo="" ):

		# The first attempt to get the histogram is only to get the minimal
		# and maximal value on the x-axis, for not predefined binning
		h_gamma = getHisto( gSignalTree, plot )
		h_fo = getHisto( foSignalTree, plot, cut="1" )
		xMin, xMax = getXMinXMax( [ h_gamma, h_fo ] )

		h_gamma = getHisto( gSignalTree, plot, color=1, firstBin=xMin,lastBin=xMax )
		h_fo = getHisto( foSignalTree, plot, cut="1", weight="weight*w_qcd", color=46, firstBin=xMin,lastBin=xMax )
		h_fo_error = getQCDErrorHisto( foSignalTree, plot, cut="1", firstBin=xMin, lastBin=xMax )
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
		SaveAs(can, "plots","qcd_afterWeighting_%s_%s"%(datasetAffix+additionalInfo, plot) )
		ROOT.SetOwnership( hPad, False )
		ROOT.SetOwnership( ratioPad, False )

def qcdClosure( fileName, opts, cuts ):
	try:
		datasetAffix = re.match(".*slim(.*)_V.*", fileName ).groups()[0]
	except:
		try:
			datasetAffix = re.match(".*slim(.*)\.root", fileName ).groups()[0]
		except:
			datasetAffix = fileName

	# Apply the cuts and divide in control and signal region
	signalCut = "met > 100"
	controlCut = "met <= 100"
	foCut = "Max$(photon.hadTowOverEm) < 0.05" \
			+"&& Max$(photon.chargedIso) < %s"%cuts["chIso"] \
			+"&& Max$(photon.neutralIso- %s *photon.pt)< %s "%(cuts["nIsoRel"], cuts["nIso"]) \
			+"&& Max$(photon.photonIso- %s *photon.pt) < %s "%(cuts["pIsoRel"], cuts["pIso"]) \
			+"&& @photon.size()>0" # +"&& ( Min$(photon.chargedIso)>2.6 || Min$(photon.sigmaIetaIeta)>0.012)"

	if True:
		foCutValentina = "Max$(photon.hadTowOverEm) < 0.05" \
			+"&& Max$(photon.chargedIso) < 15" \
			+"&& Max$(photon.neutralIso- 0.04 *photon.pt)<3.5" \
			+"&& Max$(photon.photonIso- 0.005 *photon.pt) < 1.3 " \
			+"&& ( Min$(photon.chargedIso)>2.6 || Min$(photon.sigmaIetaIeta)>0.012)" \
			+"&& @photon.size()>0"
		foCut = foCutValentina

	gTree = readTree( fileName, "photonTree" ).CopyTree("1")
	gControlTree = readTree( fileName, "photonTree").CopyTree( controlCut )
	gSignalTree = readTree( fileName, "photonTree").CopyTree( signalCut )
	foTree = readTree( fileName, "photonJetTree" ).CopyTree( foCut )
	foControlTree = readTree( fileName, "photonJetTree").CopyTree( controlCut+"&&"+foCut )
	foSignalTree = readTree( fileName, "photonJetTree").CopyTree( signalCut+"&&"+foCut )

	# Clean jets from fake object
	jets = ROOT.std.vector("tree::Jet")()
	foControlTree = clearTreeFromJets( foControlTree, jets )
	foSignalTree = clearTreeFromJets( foSignalTree, jets )
	foTree = clearTreeFromJets( foTree, jets )

	weights = getWeightHisto2D( gControlTree, foControlTree, datasetAffix )
	writeWeight2DToFile( fileName, foSignalTree, weights, "foSignalWeights" )
	foSignalTree.AddFriend( "foSignalWeights", fileName )
	writeWeight2DToFile( fileName, foTree, weights, "foWeights" )
	foTree.AddFriend( "foWeights", fileName )
	writeWeight2DToFile( fileName, foControlTree, weights, "foControlWeights" )
	foControlTree.AddFriend( "foControlWeights", fileName )

	plots = [ "met", "ht", "photon[0].ptJet", "Length$(jet.pt)" ]

	infoControl = PlotCaption()
	infoControl.controlCut()
	infoSignal = PlotCaption()
	infoSignal.signalCut()
	info = PlotCaption()
	can = ROOT.TCanvas()
	can.cd()
	for plot in plots:
		drawClosure( plot, gSignalTree, foSignalTree, can, infoSignal, datasetAffix, "_signal" )
		drawClosure( plot, gControlTree, foControlTree, can, infoControl, datasetAffix, "_control" )
		drawClosure( plot, gTree, foTree, can, info, datasetAffix )

	for plot in plots:
		drawBeforeClosure( plot, gSignalTree, foSignalTree, can, infoSignal, datasetAffix, "_signal" )
		drawBeforeClosure( plot, gControlTree, foControlTree, can, infoControl, datasetAffix, "_control" )
		drawBeforeClosure( plot, gTree, foTree, can, info, datasetAffix )


	# Delete the trees from memory, important for many iterations.
	ROOT.SetOwnership( foSignalTree, False )
	ROOT.SetOwnership( foControlTree, False )
	for tree in [ gSignalTree, gControlTree, foSignalTree, foControlTree ]:
		tree.Delete()
		del tree

def drange( start, stop, step ):
	"""Similar to 'range' function, but works with floating point steps as well."""
	r = start
	while r < stop:
		yield r
		r += step

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Calculate weighting "
			+"factors for QCD background estimation." )
	arguments.add_argument("--input", default=[""], nargs="+" )
	arguments.add_argument("--plot", default="met" )
	arguments.add_argument("--save", action="store_true", help="Save canvas as pdf.")
	opts = arguments.parse_args()

	for inName in opts.input:
		cuts = { "chIso":  4.6,
				"nIsoRel": 0.04,
				"nIso":    4.5,
				"pIsoRel": 0.005,
				"pIso":    1.3 }

		qcdClosure( inName, opts, cuts )
		for chIso in drange( 2.6, 8.2, 0.4 ):
			cuts["chIso"] = chIso
		#for nIsoRel in drange( .005, 0.024, 0.001 ):
		#	cuts["pIsoRel"] = nIsoRel
		#	for nIso in drange( 1.3, 13.5, .5 ):
		#		cuts["pIso"] = nIso

