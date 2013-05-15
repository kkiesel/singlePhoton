#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
ROOT.gSystem.Load("libTreeObjects.so")
import argparse
from multiplot import *
from treeFunctions import *
import Styles
import ratios
Styles.tdrStyle()
import os

import ConfigParser
axisConf = ConfigParser.SafeConfigParser()
axisConf.read("axis.cfg")

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="photon.pt" )
	arguments.add_argument( "--input", default="slimEWK_V01.12_tree.root" )
	arguments.add_argument( "--savePrefix", default="new" )
	opts = arguments.parse_args()

	ROOT.gROOT.SetBatch()

	slimFileName = opts.input.replace( os.path.basename(opts.input), "slim"+os.path.basename(opts.input))

	genE = Dataset( opts.input, "susyTree", "Max$(abs(genElectron.eta)) < 1.4442 && @genElectron.size() > 0", "e_{gen}", 1 )
	recE = Dataset( slimFileName, "photonElectronTree", "photon.genInformation==1", "e", 2 )
	gamma = Dataset( slimFileName, "photonTree", "photon.genInformation==1", "#gamma", 3 )

	multihisto = Multihisto()

	label, unit = readAxisConf( opts.plot, axisConf)

	histForRatio = {}

	for dataset in [genE, recE, gamma]:
		if dataset.label == "e_{gen}":
			try:
				plot = opts.plot.replace("photon","genElectron")
			except:
				plot = opts.plot
		else:
			plot = opts.plot

		hist = createHistoFromTree( dataset.tree, plot, "weight*(%s)"%(dataset.additionalCut), firstBin=10, lastBin=200 )
		hist.SetTitle(";%s%s;Entries"%(label,unit))
		hist.SetLineColor( dataset.color )
		hist.SetLineWidth(2)
		multihisto.addHisto( hist, dataset.label )
		histForRatio[dataset.label] = hist

	can = ROOT.TCanvas()

	hPad = ROOT.TPad("hPad", "Histogram", 0, 0, 1, 1)
	#hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
	hPad.cd()
	multihisto.Draw("hist")

	ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
	ratioPad.cd()
	ratioPad.SetLogy(0)
	ratioGraph = ratios.RatioGraph(histForRatio["e"], histForRatio["e_{gen}"])
	ratioGraph.draw(ROOT.gPad, yMin=0.5, yMax=1.5, adaptiveBinning=False, errors="x")
	ratioGraph.hAxis.SetYTitle( "e/e_{gen}")
	ratioGraph.graph.Draw("same p")

	can.cd()
	hPad.Draw()
	#ratioPad.Draw()
	can.SaveAs("plots/%s_%s_%s.pdf"%(opts.input[0:5],opts.plot.replace(".",""),opts.savePrefix))
	del can

