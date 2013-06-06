#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
ROOT.gSystem.Load("libTreeObjects.so")
import argparse
from multiplot import *
from treeFunctions import *
import Styles
style = Styles.tdrStyle()
style.SetOptLogy(0)
ROOT.TGaxis.SetMaxDigits(3)
import os

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

	#ratioHists[0].Divide( ratioHists[1] ) # normal division
	ratioHists[0].Divide( ratioHists[0], ratioHists[1], 1, 1, "B" ) # division using bayes theorem
	ratioHists[0].SetTitle(";%s%s;%s/%s"%(label,unit,d1.label,d2.label))
	return ratioHists[0]


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Simple EWK" )
	arguments.add_argument( "--plot", default="photon.pt" )
	arguments.add_argument( "--input", default="EWK_V01.12_tree.root" )
	arguments.add_argument( "--savePrefix", default="new" )
	opts = arguments.parse_args()

	ROOT.gROOT.SetBatch()
	import re
	# dataset name is from beginning till first '_'
	slimFileName = opts.input.replace( os.path.basename(opts.input), "slim"+os.path.basename(opts.input))
	dataset = re.match("slim([^_]*)_.*", slimFileName ).groups()[0]


	genE = Dataset( slimFileName, "genElectronTree", "1", "e_{gen}", 3 )
	genE_with_match = Dataset( slimFileName, "genElectronTree", "genElectron.phi > 4", "e_{gen, match}", 1 )
	gamma = Dataset( slimFileName, "photonTree", "1", "#gamma", 1 )
	gamma_with_match = Dataset( slimFileName, "photonTree", "photon.genInformation == 1", "#gamma_{match}", 1 )

	label, unit, binning = readAxisConf( opts.plot, axisConf )
	e_match_e_reco = divideDatasets( gamma_with_match, genE, label, unit )
	e_match = divideDatasets( genE_with_match, gamma, label, unit )


	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(0)

	datasetLabel = ROOT.TPaveText(.4,.9,.6,.98, "ndc")
	datasetLabel.SetFillColor(0)
	datasetLabel.SetBorderSize(0)
	datasetLabel.AddText( dataset )


	e_match_e_reco.GetYaxis().SetTitle("#varepsilon_{match}#upointf_{e_{gen}#rightarrow#gamma}")
	e_match_e_reco.Draw("e")
	datasetLabel.Draw()
	can.SaveAs("plots/%sEfficiencyFakeRate.pdf"%dataset)
	e_match.GetYaxis().SetTitle("#varepsilon_{match}")
	e_match.Draw("e")
	datasetLabel.Draw()
	can.SaveAs("plots/%sEfficiency.pdf"%dataset)

	h = e_match_e_reco.Clone("fakerate")
	h.GetYaxis().SetTitle("f_{e_{gen}#rightarrow #gamma}")

	h.Divide( e_match )
	#h = divideDatasets( gamma, recE, label, unit )


	yutaro = h.Clone("yutaro")
	yutaro.SetBinContent(1,0.0131)
	yutaro.SetBinError  (1,0.0004)
	yutaro.SetBinContent(2,0.0146)
	yutaro.SetBinError  (2,0.0002)
	yutaro.SetBinContent(3,0.0148)
	yutaro.SetBinError  (3,0.0001)
	yutaro.SetBinContent(4,0.0111)
	yutaro.SetBinError  (4,0.0002)
	yutaro.SetBinContent(5,0.0111)
	yutaro.SetBinError  (5,0.0004)
	yutaro.SetBinContent(6,0.0085)
	yutaro.SetBinError  (6,0.0005)
	yutaro.SetLineColor(2)
	yutaro.SetMarkerColor(2)

	h.SetMaximum(max(h.GetMaximum(),yutaro.GetMaximum())+0.002)
	h.SetMinimum(min(h.GetMinimum(),yutaro.GetMinimum())-0.002)
	h.Draw("e")
	yutaro.Draw("same e0")

	leg = myLegend(.5,.70,.95,.92)
	leg.AddEntry( h, h.GetYaxis().GetTitle(), "lp")
	leg.AddEntry( yutaro, "Yutaro's f_{e#rightarrow#gamma}", "lp" )
	leg.SetBorderSize(1)
	leg.Draw()

	datasetLabel.Draw()


	saveName = "%s_%s_%s_%s"%(h.GetYaxis().GetTitle(),dataset,opts.plot,opts.savePrefix)
	saveName = saveName.replace("/","VS")
	saveName = saveName.replace(" ","_")
	unallowedCharacters = ["{","}","(",")","#","|","."]
	for char in unallowedCharacters:
		saveName = saveName.replace( char, "" )

	can.SaveAs("plots/%s.pdf"%saveName)


