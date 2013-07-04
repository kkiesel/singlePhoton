#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
import argparse
import ConfigParser
from math import sqrt
from multiplot import *
from treeFunctions import *
import ratios
import Styles

Styles.tdrStyle()
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libTreeObjects.so")

# to use user defined help message, sys.arv has to be sent to python and not
# to TApplication
ROOT.PyConfig.IgnoreCommandLineOptions = True

def writeWeight1D( fileName, tree, h_weight, variable ):
	f = ROOT.TFile( fileName, "update" )

	weightTree = ROOT.TTree("weightTree", "Tree containing QCD weights" )
	import numpy
	weight = numpy.zeros( 1, dtype=float)
	weight_error = numpy.zeros( 1, dtype=float)

	weightTree.Branch( "w_qcd", weight, "w_qcd/D" )
	weightTree.Branch( "w_qcd_error", weight_error, "w_qcd_error/D" )

	variable = variable.replace("[0]",".at(0)")
	variable = "tree."+variable

	for event in tree:
		if not event.GetReadEntry()%10000:
			print "%s / %s"%(event.GetReadEntry(), event.GetEntries() )
		pt = eval( variable )
		weight[0] = h_weight.GetBinContent( h_weight.FindBin( pt ) )
		weight_error[0] = h_weight.GetBinError( h_weight.FindBin( pt ) )
		weightTree.Fill()

	f.cd()
	weightTree.Write()
	f.Close()

def writeWeight2D( fileName, tree, h_weight, weightTreeName ):
	weightTree = ROOT.TTree( weightTreeName, "Tree containing QCD weights" )
	import numpy
	weight = numpy.zeros( 1, dtype=float)
	weight_error = numpy.zeros( 1, dtype=float)

	weightTree.Branch( "w_qcd", weight, "w_qcd/D" )
	weightTree.Branch( "w_qcd_error", weight_error, "w_qcd_error/D" )

	for event in tree:
		if not event.GetReadEntry()%10000:
			print "%s / %s"%(event.GetReadEntry(), event.GetEntries() )
		b = h_weight.FindBin( event.photon.at(0).ptJet, event.jet.size() )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()

	f = ROOT.TFile( fileName, "update" )
	f.cd()
	weightTree.Write()
	f.Close()

def drawDifferentWeights1D( weight2D, cutText ):
	leg = ROOT.TLegend(.2, .6, .5, .9 )
	leg.SetHeader("n_{Jet}")
	weightHistos = []
	for yBin in range( weight2D.GetNbinsY()+1 ):
		h = weight2D.ProjectionX( "_px%s"%yBin, yBin, yBin+1 )
		h.SetLineColor( yBin+1 )
		h.SetMarkerColor( h.GetLineColor() )
		h.SetTitle("w")
		weightHistos.append( h )

	can = ROOT.TCanvas()
	can.SetLogy(0)
	can.cd()

	weightHistos[0].Draw()
	for h in weightHistos[1:]:
		h.Draw("same")

	cutText.Draw()
	can.SaveAs("qcd_preWeight_weights1D.pdf")
	return

def writeWeights( fileName, gControlTree, foControlTree, foSignalTree, plotStuff=None ):
	cutText, datasetAffix, cutSaveAffix = plotStuff

	xlabel, xunit, xbinning = readAxisConf( "photon[0].ptJet" )
	ylabel, yunit, ybinning = "n_{Jet}", "", [ 1.5, 2.5 ]

	weight_numerator = createHistoFromTree2D( gControlTree, "Length$(jet.pt):photon[0].ptJet", "weight", xbinning, ybinning )
	weight_denominator = createHistoFromTree2D( foControlTree, "Length$(jet.pt):photon[0].ptJet", "weight", xbinning, ybinning )
	weight2D = divideHistos( weight_numerator, weight_denominator )

	drawDifferentWeights1D( weight2D, cutText )

	Styles.tdrStyle2D()
	can2D = ROOT.TCanvas()
	can2D.cd()
	for hist, name in [ (weight_numerator,"numerator"), (weight_denominator, "denominator"), (weight2D, "weight2D") ]:
		h = appendFlowBin2D( hist, 10, 10, 0, .53 )
		h.Draw("colz")
		cutText.Draw()
		can2D.SaveAs("qcd_preWeight_%s_%s_%s.pdf"%(name,datasetAffix,cutSaveAffix))
	Styles.tdrStyle()

	writeWeight2D( fileName, foControlTree, weight2D, "foControlWeights" )
	writeWeight2D( fileName, foSignalTree, weight2D, "foSignalWeights" )

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

def qcdClosure( fileName, opts, cuts ):
	# process cut string for readability
	cutSaveAffix = "%s_%s_%s_%s_%s"%(cuts["chIso"], cuts["nIsoRel"], cuts["nIso"], cuts["pIsoRel"], cuts["pIso"])
	cutText = ROOT.TPaveText(0,.95, 1, .99, "ndc" )
	cutText.AddText( "#sigma_{i#etai#eta}<0.014  Iso^{#pm}<%s  Iso^{0}<%s+%sp_{T}  Iso^{#gamma}<%s+%sp_{T}"%(cuts["chIso"],cuts["nIso"],cuts["nIsoRel"],cuts["pIso"],cuts["pIsoRel"]) )
	cutText.SetFillColor(0)
	cutText.SetLineWidth(0)
	cutText.SetLineColor(0)
	cutText.SetShadowColor(0)
	import re
	datasetAffix = re.match(".*slim(.*)_V.*", fileName ).groups()[0]

	# Apply the cuts and divide in control and signal region
	signalCut = "met > 100"
	controlCut = "met <= 100"
	foCut = "Min$(photon.pt) > 80" \
			+"&& Max$(photon.sigmaIetaIeta) < 0.014" \
			+"&& Max$(photon.hadTowOverEm) < 0.05" \
			+"&& Max$(photon.chargedIso) < 15" \
			+"&& Max$(photon.neutralIso-0.04*photon.pt)<3.5" \
			+"&& Max$(photon.photonIso-0.005*photon.pt) < 1.3" \
			+"&& ( Min$(photon.chargedIso)>2.6 || Min$(photon.sigmaIetaIeta)>0.012 )" \
			+"&& @photon.size()>0"

	foCut = "Max$(photon.hadTowOverEm) < 0.05" \
			+"&& Max$(photon.chargedIso) < %s"%cuts["chIso"] \
			+"&& Max$(photon.neutralIso- %s *photon.pt)< %s "%(cuts["nIsoRel"], cuts["nIso"]) \
			+"&& Max$(photon.photonIso- %s *photon.pt) < %s "%(cuts["pIsoRel"], cuts["pIso"]) \
			+"&& @photon.size()>0" # +"&& ( Min$(photon.chargedIso)>2.6 || Min$(photon.sigmaIetaIeta)>0.012)"

	gSignalTree = readTree( fileName, "photonTree").CopyTree( signalCut )
	gControlTree = readTree( fileName, "photonTree").CopyTree( controlCut )
	foSignalTree = readTree( fileName, "photonJetTree").CopyTree( signalCut+"&&"+foCut )
	foControlTree = readTree( fileName, "photonJetTree").CopyTree( controlCut+"&&"+foCut )

	# Clean jets from fake object
	jets = ROOT.std.vector("tree::Jet")()
	foControlTree = clearTreeFromJets( foControlTree, jets )
	foSignalTree = clearTreeFromJets( foSignalTree, jets )

	writeWeights( fileName, gControlTree, foControlTree, foSignalTree, [cutText, datasetAffix, cutSaveAffix] )

	foSignalTree.AddFriend( "foSignalWeights", fileName )
	foControlTree.AddFriend( "foControlWeights", fileName )

	plots = [ "met", "ht", "photon[0].ptJet", "Length$(jet.pt)" ]

	can = ROOT.TCanvas()
	can.cd()
	for plot in plots:
		# The first attempt to get the histogram is only to get the minimal
		# and maximal value on the x-axis, for not predefined binning
		h_gamma = getHisto( gSignalTree, plot )
		h_fo = getHisto( foSignalTree, plot, cut="1" )
		xMin, xMax = getXMinXMax( [ h_gamma, h_fo ] )

		# for integers, adjust nBins and shift by 0.5
		if "Length$(" in plot or "nVertex" == plot:
			xMin -= .5
			xMax += .5
			nBins = int(xMax-xMin)
		else:
			nBins = 20

		# there is a muon with 7000 GeV, which destroys the automatic binning
		if plot == "muon.pt":
			xMax = 200

		h_gamma = getHisto( gSignalTree, plot, color=1, nBins=nBins, firstBin=xMin,lastBin=xMax )
		h_fo = getHisto( foSignalTree, plot, cut="1", weight="weight*w_qcd", color=46, nBins=nBins, firstBin=xMin,lastBin=xMax )
		h_fo_error = getQCDErrorHisto( foSignalTree, plot, cut="1",nBins=nBins, firstBin=xMin, lastBin=xMax )
		h_fo_error.SetFillColor(2)
		h_fo_error.SetFillStyle(3254)
		h_fo_error.SetMarkerSize(0)

		muhisto = Multihisto()
		muhisto.leg.SetHeader( datasetAffix )
		muhisto.addHisto( h_gamma, "#gamma", draw="hist e0" )
		muhisto.addHisto( h_fo, "#gamma_{jet}#upointw", draw="hist e0")
		muhisto.addHisto( h_fo_error, "#sigma_{w}", draw="e2")

		hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
		hPad.cd()
		muhisto.Draw()
		cutText.Draw()

		ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
		ratioPad.cd()
		ratioPad.SetLogy(0)
		ratioGraph = ratios.RatioGraph(h_fo, h_gamma)
		ratioGraph.draw(ROOT.gPad, yMin=0.5, yMax=1.5, adaptiveBinning=False, errors="yx")
		ratioGraph.graph.Draw("same p e0") # draw nice points
		ratioGraph.hAxis.SetYTitle( "#gamma_{pred}/#gamma")
		can.cd()
		hPad.Draw()
		ratioPad.Draw()
		SaveAs(can, "plots","qcd_%s_%s_afterWeighting_%s"%(datasetAffix, plot,cutSaveAffix))
		ROOT.SetOwnership( hPad, False )
		ROOT.SetOwnership( ratioPad, False )

	# Delete the trees from memory, important for many iterations.
	for tree in [ gSignalTree, gControlTree ]:
		tree.Delete()
		del tree
	del foSignalTree
	del foControlTree


def drange( start, stop, step ):
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
		cuts = { "chIso":  2.6,
				"nIsoRel": 0.04,
				"nIso":    3.5,
				"pIsoRel": 0.005,
				"pIso":    1.3 }

		for chIso in drange( 2.6, 7.2, 0.4 ):
			cuts["chIso"] = chIso
			qcdClosure( inName, opts, cuts )
		#for nIsoRel in drange( .005, 0.024, 0.001 ):
		#	cuts["pIsoRel"] = nIsoRel
		#	for nIso in drange( 1.3, 13.5, .5 ):
		#		cuts["pIso"] = nIso

