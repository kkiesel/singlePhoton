#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from predictions import *
from qcdClosure import drawWeightHisto

ROOT.gStyle.SetOptLogy(0)

def getCombinatoricalBkg( filenames, plot ):
	# calculates invariant masses between all muons (also from different events)

	tree = ROOT.TChain("photonTree")
	for filename in filenames[0:1]:
		tree.AddFile( filename )

	# 1st loop to get all muons
	import copy
	muons = []
	for event in tree:
		for muon in event.muons:
			muons.append( copy.copy(muon) )

	# create empty histogram
	histo = getHisto( tree, plot, cut="0" )
	from calculateAdditionalVariables import M
	for event in tree:
		if event.muons.size() == 1:
			for muon in muons:
				minv = M( event.muons.at(0), muon )
				histo.Fill( minv )

	return histo


def drawChi2( tuples, plot ):

	# draw graph
	gr = ROOT.TGraph()
	gr.SetTitle(";k-factor;#chi^{2}/ndf")
	for i, (k, chi2) in enumerate( tuples ):
		gr.SetPoint(i, k, chi2 )
	gr.Draw("ap")

	#getMinima
	minkFactor, minChi2 = min( tuples, key=lambda t: t[1] )

	# draw lines
	l = ROOT.TLine( gr.GetXaxis().GetXmin(), minChi2+1, gr.GetXaxis().GetXmax(), minChi2+1 )
	l.Draw()

	# get interception with graph
	interceptions = []
	for i in range(len(tuples)-1):
		if minChi2+1 >= tuples[i][1] and minChi2+1 < tuples[i+1][1] \
			or minChi2+1 <= tuples[i][1] and minChi2+1 > tuples[i+1][1]:
			interceptions.append( tuples[i] )

	intLine = ROOT.TLine()
	intLine.SetLineStyle(7)
	for interception in interceptions:
		intLine.DrawLine( interception[0], gr.GetYaxis().GetXmin(), interception[0], minChi2+1 )

	text = "k-factor = %s"%minkFactor
	if len(interceptions) > 0:
		text += "_{-%s}"%(minkFactor-interceptions[0][0])
	if len(interceptions) > 1:
		text += "^{+%s}"%(interceptions[1][0]-minkFactor)


	# draw text
	kText = ROOT.TLatex(0.4, .8, text )
	kText.SetNDC()
	kText.Draw()
	if len(interceptions) and False:
		up = minkFactor-interceptions[0][0]
		down = interceptions[1][0]-minkFactor if len(interceptions)>1 else up
		relUncert = (up+down)/(2*minkFactor)*100

		skText = ROOT.TLatex(0.5, .7, "rel uncert = %i%%"%relUncert )
		skText.SetNDC()
		skText.Draw()

	ROOT.gPad.SaveAs("plots/chi2Minimization_%s.pdf"%plot)


def getkFactor( dataFiles, bkgFiles, plot, cut ):

	kFactorMin = 0.5
	kFactorMax = 4
	kFactorN = 100

	data = getHists( dataFiles, plot, cut )
	bkg = getHists( bkgFiles, plot, cut )

	tuples = []

	for i in range( kFactorN ):
		kFactor = kFactorMin + 1.*(kFactorMax-kFactorMin)/kFactorN*i
		scaledbkg = bkg.Clone(randomName())
		scaledbkg.Scale( kFactor )
		data.Draw()
		scaledbkg.Draw("hist same")
		chi2 = 0
		ndf = 0
		for bin in range( data.GetNbinsX()+2 ):
			chi2 += (data.GetBinContent(bin)-scaledbkg.GetBinContent(bin))**2 / (data.GetBinError(bin)**2 + scaledbkg.GetBinError(bin)**2 ) if data.GetBinError(bin) or scaledbkg.GetBinError(bin) else 0
			if data.GetBinContent(bin) and scaledbkg.GetBinContent(bin):
				ndf += 1
		if ndf: chi2 /= ndf

		#chi2 = data.Chi2Test( scaledbkg, "uw CHI2/NDF p" )
		#ROOT.gPad.SaveAs("test%s.pdf"%i)
		#print kFactor, chi2
		tuples.append( (kFactor, chi2) )

	drawChi2( tuples, plot )

	return min( tuples, key=lambda t: t[1] )[0]


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("--plot", default = "mmm" )
	opts = arguments.parse_args()

	cut = "std::abs(photons[0].eta)<1.4442 && met<100"
	# cleaning:
	cut += " && {} > 10".format(opts.plot)
	chi2Cut = cut + "&& 60<{0} && {0}<120".format(opts.plot)

	treeVersion = 24
	dataFiles = [ "PhotonHad%s_V03.%s_tree.root"%(x,treeVersion) for x in ["A","B","C","D" ] ]
	data = getHists( dataFiles, opts.plot, cut )

	bkg = {}
	bkg["gjets"] = { "files": ["slimGJets_400_inf_V03.%s_tree.root"%treeVersion, "slimGJets_200_400_V03.%s_tree.root"%treeVersion ], "title":"#gammaJet", "color": ROOT.kCyan }

	bkg["zll"] = { "files": ["slimZGammaLL_V02.%s_tree.root"%treeVersion], "title": "#gammaZ#rightarrow#gammall", "color": 2 }
	bkg["tt"] = { "files": ["slimTTGamma_V03.%s_tree.root"%treeVersion], "title": "#gammat#bar{t}", "color": ROOT.kBlue }
	bkg["wjets"] = { "files": ["slimWJets_250_300_V03.24_tree.root", "slimWJets_300_400_V03.24_tree.root", "slimWJets_400_inf_V03.24_tree.root" ], "title": "W", "color": ROOT.kGreen+4 }
	bkg["wgamma"] = { "files": ["slimWGamma_130_inf_V03.24_tree.root", "slimWGamma_50_130_V03.24_tree.root" ], "title": "#gammaW", "color": ROOT.kGreen-4 }
	#bkg["qcd"] = { "files": ["slimQCD_1000_inf_V03.24_tree.root", "slimQCD_250_500_V03.24_tree.root", "slimQCD_500_1000_V03.24_tree.root"], "title":"Multijet", "color": ROOT.kCyan+3 }

	nestedBkgFiles = [ bkg[a]["files"] for a in bkg.keys()]
	bkgFiles = [item for sublist in nestedBkgFiles for item in sublist]

	kFactor = getkFactor( dataFiles, bkgFiles, opts.plot, chi2Cut )

	#signal = getHists( ["slimW_1700_720_375_V03.24_tree.root" ], opts.plot, cut )
	#signal.SetLineColor( ROOT.kGreen )
	#signal.SetLineWidth(2)

	mh = Multihisto()
	mh.setMinimum(0)
	mh.addHisto( data, "Data", draw="pe" )
	for name, d in bkg.iteritems():
		histo = getHists( d["files"], opts.plot, cut )
		histo.SetLineColor( d["color"] )
		histo.Scale( kFactor )
		mh.addHisto( histo, d["title"], True )

	#mh.addHisto( signal, "Wino", False )

	#combiBkg = getCombinatoricalBkg( dataFiles, opts.plot )
	#combiBkg.SetLineWidth(2)
	#combiBkg.SetLineColor( ROOT.kBlue )
	#combiBkg.Scale( data.Integral(0, data.FindBin(70), "width") / combiBkg.Integral(0,data.FindBin(70),"width"))
	#mh.addHisto( combiBkg, "bkg", draw="hist e" )


	mh.Draw()

	ROOT.gPad.SaveAs( "plots/isrkFactor_%s.pdf"%opts.plot )

