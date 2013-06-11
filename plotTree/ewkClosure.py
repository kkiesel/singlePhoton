#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
import re
import argparse
from math import sqrt
from multiplot import *
from treeFunctions import *
import ratios
import Styles

Styles.tdrStyle()
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libTreeObjects.so")

def applyFakeRate( histo, f, e_f ):
	for i in range( histo.GetNbinsX() +1 ):
		binContent = histo.GetBinContent(i)
		histo.SetBinContent( i, histo.GetBinContent(i)*f)
		# sigma_{ef} = sqrt( ( e*e_f )**2 + ( e_e*f )**2 )
		histo.SetBinError( i, sqrt((binContent*e_f)**2 + (histo.GetBinError(i)*f)**2 ) )
	return histo

def closure( fileName, opts ):
	##################################
	fakeRate = 0.0084
	fakeRateError = 0.0006 # stat

	# correct fake rate, if it is estimated with yutaros method
	fakeRateError = fakeRateError / (1-fakeRate)**2
	fakeRate = fakeRate / ( 1 - fakeRate )

	##################################

	# dataset name is from beginning till first '_'
	datasetAffix = re.match(".*slim([^_]*)_.*", fileName ).groups()[0]

	label, unit, binning = readAxisConf( opts.plot )

	recE = extractHisto( Dataset( fileName, "photonElectronTree", "Max$(photon.isGenElectron()) && Min$(photon.pt)>80", "e", 2 ), opts.plot)
	recE.SetMarkerSize(0)
	recE.SetFillColor( recE.GetLineColor() )
	recE.SetFillStyle(3254)

	recE = applyFakeRate( recE, fakeRate, fakeRateError )

	gamma = extractHisto( Dataset( fileName, "photonTree", "Max$(photon.isGenElectron()) && Min$(photon.pt)>80", "#gamma", 1 ), opts.plot)
	gamma.SetMarkerSize(0)

	multihisto = Multihisto()
	multihisto.leg.SetHeader( datasetAffix )
	multihisto.addHisto( recE, "e#upoint#tildef_{e#rightarrow#gamma}", draw="e2" )
	multihisto.addHisto( gamma, "#gamma", draw="e hist" )

	can = ROOT.TCanvas()

	hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
	hPad.cd()
	multihisto.Draw()

	ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
	ratioPad.cd()
	ratioPad.SetLogy(0)
	ratioGraph = ratios.RatioGraph(gamma, recE)
	ratioGraph.draw(ROOT.gPad, yMin=0.5, yMax=1.5, adaptiveBinning=False, errors="yx")
	ratioGraph.hAxis.SetYTitle( "#gamma/(e#upoint#tilde{f})")

	"""
	scale_for_points = 3.66 # mean weight for w-jets photonTree
	for point in range( ratioGraph.graph.GetN() ):
		if not ratioGraph.graph.GetErrorY(point):
			ratioGraph.graph.SetPointEYhigh( point, 1.14*scale_for_points/ ratioGraph.denominator.GetBinContent(point+1) )
	"""

	ratioGraph.graph.Draw("same p")

	# draw systematic uncertanty in ratio:
	syst = recE.Clone( "hist_ratio_syst" )
	for bin in range( syst.GetNbinsX()+1 ):
		if syst.GetBinContent(bin):
			syst.SetBinError( bin, syst.GetBinError(bin) / syst.GetBinContent(bin) )
		syst.SetBinContent( bin, 1 )
	syst.Draw("same e2")

	can.cd()
	hPad.Draw()
	ratioPad.Draw()
	can.SaveAs("plots/%s_%s_%s.pdf"%(datasetAffix,opts.plot.replace(".",""),opts.savePrefix))
	del can


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--input", default=["EWK_V01.12_tree.root"], nargs="+" )
	arguments.add_argument( "--savePrefix", default="new" )
	opts = arguments.parse_args()

	for inName in opts.input:
		closure( inName, opts )

