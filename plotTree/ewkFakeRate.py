#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
ROOT.gSystem.Load("libTreeObjects.so")
import argparse
from multiplot import *
from treeFunctions import *
import Styles
Styles.tdrStyle()

import ConfigParser
axisConf = ConfigParser.SafeConfigParser()
axisConf.read("axis.cfg")

def divideDatasets( d1, d2, label, unit ):
	# d1, d2 are multiplot.Datasets (own definition)
	ratioHists = []
	for dataset in [d1, d2]:
		if dataset.label == "e_{gen}":
			try:
				plot = opts.plot.replace("photon","genElectron")
			except:
				plot = opts.plot
		else:
			plot = opts.plot

		yutarosBinning = [ 25, 35, 40, 50, 60, 80, 100 ]
		hist = createHistoFromTree( dataset.tree, plot, "weight*(%s)"%(dataset.additionalCut), nBins=yutarosBinning)
		hist.SetTitle(";%s%s;Entries"%(label,unit))
		hist.SetLineColor( dataset.color )
		hist.SetLineWidth(2)
		ratioHists.append( hist )

	ratioHists[0].Divide( ratioHists[1] )
	ratioHists[0].SetTitle(";%s%s;%s/%s"%(label,unit,d1.label,d2.label))
	return ratioHists[0]


if True:
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="photon.pt" )
	arguments.add_argument( "--dataset", default="EWK" )
	opts = arguments.parse_args()

	version = "07"
	fileName = "slim%s_V01.%s_tree.root"%(opts.dataset, version )

	genE = Dataset( fileName, "genElectronTree", "abs(genElectron.eta) < 1.479", "e_{gen}", 3 )
	recE = Dataset( fileName, "photonElectronTree", "@genElectron.size()>0 && photon[0].pixelseed < 0", "e", 2 )
	gamma = Dataset( fileName, "photonTree", "@genElectron.size()>0 && photon[0].pixelseed < 0", "#gamma", 1 )

	label, unit = readAxisConf( opts.plot, axisConf )
	h = divideDatasets( gamma, recE, label, unit )


	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(0)

	yutaro = h.Clone("yutaro")
	yutaro.SetBinContent(1,0.0196)
	yutaro.SetBinError(1,0.0004)
	yutaro.SetBinContent(2,0.0198)
	yutaro.SetBinError(2,0.0006)
	yutaro.SetBinContent(3,0.0202)
	yutaro.SetBinError(3,0.0002)
	yutaro.SetBinContent(4,0.0180)
	yutaro.SetBinError(4,0.0004)
	yutaro.SetBinContent(5,0.0162)
	yutaro.SetBinError(5,0.0006)
	yutaro.SetBinContent(6,0.0148)
	yutaro.SetBinError(6,0.0010)
	yutaro.SetLineColor(2)
	yutaro.SetMarkerColor(2)

	h.SetMaximum(max(h.GetMaximum(),yutaro.GetMaximum())+0.002)
	h.SetMinimum(min(h.GetMinimum(),yutaro.GetMinimum())-0.002)
	h.Draw("e")
	yutaro.Draw("same e0")

	can.SaveAs("pt_fake_new.pdf")
