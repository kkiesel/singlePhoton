#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
ROOT.gSystem.Load("libTreeObjects.so")
import argparse
from multiplot import *
from treeFunctions import *
import Styles
from math import sqrt
import ratios
Styles.tdrStyle()

import ConfigParser
axisConf = ConfigParser.SafeConfigParser()
axisConf.read("axis.cfg")

def applyFakeRate( histo, f, e_f ):
	for i in range( histo.GetNbinsX() +1 ):
		binContent = histo.GetBinContent(i)
		histo.SetBinContent( i, histo.GetBinContent(i)*f)
		# sigma_{ef} = sqrt( ( e*e_f )**2 + ( e_e*f )**2 )
		histo.SetBinError( i, sqrt((binContent*e_f)**2 + (histo.GetBinError(i)*f)**2 ) )
	return histo

def extractHisto( dataset, plot ):
	histo = createHistoFromTree( dataset.tree, plot, "weight*(%s)"%(dataset.additionalCut), nBins=30, firstBin=0, lastBin=200)
	label, unit = readAxisConf( opts.plot, axisConf )
	if unit:
		histo.SetTitle(";%s%s;Entries / %s"%(label, unit, unit[2:-1]))
	else:
		histo.SetTitle(";%s%s;Entries"%(label, unit))
	histo.SetLineColor( dataset.color )
	histo.SetMarkerColor( dataset.color )
	histo.SetMarkerSize(0)
	return histo

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--input", default="EWK_V01.12_tree.root" )
	arguments.add_argument( "--savePrefix", default="new" )
	opts = arguments.parse_args()


	##################################
	fakeRate = 0.0084 # dy->ee mc yutaro
	fakeRateError = 0.0005 # stat

	# correct fake rate, if it is estimated with yutaros method
	fakeRateError = fakeRateError / (1-fakeRate)**2
	fakeRate = fakeRate / ( 1 - fakeRate )

	##################################

	ROOT.gROOT.SetBatch()
	import re
	# dataset name is from beginning till first '_'
	slimFileName = opts.input.replace( os.path.basename(opts.input), "slim"+os.path.basename(opts.input))
	datasetAffix = re.match("slim([^_]*)_.*", slimFileName ).groups()[0]

	recE = extractHisto( Dataset( slimFileName, "photonElectronTree", "1", "e", 2 ), opts.plot )
	recE.SetFillColor( recE.GetLineColor() )
	recE.SetFillStyle(3254)

	recE = applyFakeRate( recE, fakeRate, fakeRateError )

	gamma = extractHisto( Dataset( slimFileName, "photonTree", "Max$(photon.isGenElectron())", "#gamma", 1 ), opts.plot )

	label, unit = readAxisConf( opts.plot, axisConf)
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
	ratioGraph.graph.Draw("same p")

	# draw systematic uncertanty in ratio:
	syst = recE.Clone( "newname" )
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

