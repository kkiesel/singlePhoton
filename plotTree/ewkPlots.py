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

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--dataset", default="EWK" )
	opts = arguments.parse_args()

	version = "07"
	fileName = "slim%s_V01.%s_tree.root"%(opts.dataset, version )

	genE = Dataset( fileName, "genElectronTree", "abs(genElectron.eta) < 1.479", "e_{gen}", 1 )
	recE = Dataset( fileName, "photonElectronTree", "@genElectron.size()>0 && photon[0].pixelseed < 0", "e", 2 )
	gamma = Dataset( fileName, "photonTree", "@genElectron.size()>0 && photon[0].pixelseed < 0", "#gamma", 3 )

	multihisto = Multihisto()

	label, unit = readAxisConf( opts.plot, axisConf)

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

	can = ROOT.TCanvas()
	can.cd()
	multihisto.Draw("hist")
	can.SaveAs("pt_%s_new.pdf"%opts.dataset)


