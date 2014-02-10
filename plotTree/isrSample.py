#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

colors = {
		"GJets": ROOT.kCyan-7,
		"QCD": ROOT.kCyan+3,
		"GJets_200_400": ROOT.kCyan-10,
		"GJets_400_inf":ROOT.kCyan-7,
		"QCD_250_500": ROOT.kCyan+3,
		"QCD_500_1000": ROOT.kCyan+2,
		"QCD_1000_inf": ROOT.kCyan+1,
		"WJets": ROOT.kGreen,
		"WJets_250_300": ROOT.kGreen-1,
		"WJets_300_400": ROOT.kGreen-2,
		"WJets_400_inf": ROOT.kGreen-3,
		"TTJets": ROOT.kGreen-3,
		"WGamma": ROOT.kRed-7,
		"WGamma_50_130": ROOT.kRed-6,
		"WGamma_130_inf": ROOT.kRed-5,
		"ZGamma": ROOT.kBlue,
		"ZGammaNuNu": ROOT.kMagenta-3,
		"ZGammaNuNu1_400_inf": ROOT.kMagenta-2,
		"ZGammaNuNu2_400_inf": ROOT.kMagenta-4,
		"ZGammaNuNu3_400_inf": ROOT.kMagenta-6,
		"ZGammaLL": ROOT.kMagenta-2,
		"TTGamma": ROOT.kRed,
		"TTHadronic": ROOT.kGreen-6,
		"TTFull": ROOT.kGreen,
		"TTSemi": ROOT.kGreen-2
		}

info = PlotCaption()
info = ROOT.TLatex(0,.97, "#text{CMS Private Work }#hspace{5cm}#SI{19.8}{fb^{-1}}#, #sqrt{s}=#SI{8}{TeV}#, #geq1#gamma,#geq2#text{jets}" )
info.SetNDC()
info.SetTextSize(1)

ROOT.gStyle.SetCanvasDefH(1000)
ROOT.gStyle.SetCanvasDefW(2000)
ROOT.gStyle.SetPaperSize(14.65,50)

# Margins:
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.10)
ROOT.gStyle.SetPadRightMargin(0.02)



def getHists( filenames, plot="met", cut="1" ):
	endHist = None
	for filename in filenames:
		tree = readTree( filename, "photonTree" )
		hist = getHisto( tree, plot, color=1, fillEmptyBins=True, cut=cut )

		if endHist: endHist.Add( hist )
		else: endHist = hist

	return endHist

def compareBinnedSamples( histList1, histList2, plot="met" ):

	cut = "!@electrons.size() && !@muons.size()"
	treeName = "photonTree"

	mh = Multihisto()

	#if "ZGamma" in histList1[0]:
	#	mh.setMinimum(0.01)

	# sum all histos in histList1
	hist1 = None
	for fileName in histList1:
		datasetAbbr1 = getDatasetAbbr( fileName )
		tree = readTree( fileName, treeName )
		h = getHisto( tree, plot, cut=cut, color=colors[datasetAbbr1] )
		if hist1: hist1.Add( h )
		else: hist1 = h
	mh.addHisto( hist1, datasetToLatex(datasetAbbr1), draw="hist e" )

	# stack histos in histList2
	for fileName in histList2:
		datasetAbbr = getDatasetAbbr( fileName )
		tree = readTree( fileName, treeName )
		h = getHisto( tree, plot, cut=cut, color=colors[datasetAbbr] )
		mh.addHisto( h, datasetToLatex(datasetAbbr), toStack=True )

	mh.Draw()

	stack = mh.stack.GetStack().Last()
	stack.SetMarkerSize(0)
	stack.Draw("same e0")
	errorHist = mh.stack.GetStack().Last().Clone(randomName() )
	errorHist.SetMarkerColor(1)
	errorHist.Draw("same e")
	if plot != "met":
		datasetAbbr1 += plot
	SavePad("compareBinned"+datasetAbbr1 )


def inclusiveAndIsrSamples( fList1, fList2, saveAffix="" ):

	cut = "!@electrons.size() && !@muons.size()"
	treeName = "photonTree"
	plot = "met"

	if saveAffix == "pt130":
		cut += " && photons[0].pt>130"
		saveAffix = "_"+saveAffix

	mh = Multihisto()

	h1 = getHists( fList1, plot, cut )
	h1gen = getHists( fList1, plot, cut+"&&photons[0].isGen(0)" )

	h2 = getHists( fList2, plot, cut)
	h2gen = getHists( fList2, plot, cut+"&&photons[0].isGen(0)" )

	for h in h2, h2gen:
		h.SetLineColor(2)
	for h in h1gen, h2gen:
		h.SetLineStyle(2)

	for h in [h2,h2gen,h1,h1gen]:
		h.SetLabelSize( 1./31.5562/0.502113, "xyz" )
		h.SetTitleSize( 1./31.5562/0.502113, "xyz" )
		h.GetXaxis().SetTitleOffset(1.1)
		h.GetYaxis().SetTitleOffset(0.85)
		if plot == "met":
			h.GetXaxis().SetTitle("#met#text{ [GeV]}")

	abbr1 = shortName( fList1 )
	abbr2 = shortName( fList2 )

	mh = Multihisto()
	mh.leg.SetX1(0.6)
	mh.leg.SetY1(0.6)
	mh.addHisto( h1, datasetToLatex( abbr1 ) )
	mh.addHisto( h1gen, "#text{match to gen }#gamma" )
	mh.addHisto( h2, datasetToLatex( abbr2 ) )
	mh.addHisto( h2gen, "#text{match to gen }#gamma" )

	mh.Draw()
	info.SetTextSize(1./31.5562/0.502113)
	info.Draw()
	SavePad( "inclusiveAndIsrSample_%s"%abbr1 )

	ROOT.gPad.SaveAs("/home/knut/master/documents/thesis/plots/inclusiveAndIsrSample_%s%s.tex"%(abbr1,saveAffix) )
	correctTiksPlot( "/home/knut/master/documents/thesis/plots/inclusiveAndIsrSample_%s%s.tex"%(abbr1,saveAffix) )



if __name__ == "__main__":

	compareBinnedSamples(
		["slimWJets_V02.43_tree.root"],
		[ "slimWJets_250_300_V02.43_tree.root", "slimWJets_300_400_V02.43_tree.root", "slimWJets_400_inf_V02.43_tree.root" ]
		)

	compareBinnedSamples(
		[ "slimWGamma_V02.19_tree.root" ],
		[ "slimWGamma_50_130_V02.43_tree.root", "slimWGamma_130_inf_V02.43_tree.root" ]
		)

	compareBinnedSamples(
		[ "slimZGamma_V02.43_tree.root" ],
		[ "slimZGammaNuNu_V02.43_tree.root" ]
		#[ "slimZGammaLL_V02.43_tree.root", "slimZGammaNuNu_V02.43_tree.root" ]
		)

	compareBinnedSamples(
		[ "slimZGamma_V02.43_tree.root" ],
		[ "slimZGammaNuNu_V02.43_tree.root" ],
		#[ "slimZGammaLL_V02.43_tree.root", "slimZGammaNuNu_V02.43_tree.root" ],
		"photons[0].ptJet()"
		)

	compareBinnedSamples(
		[ "slimTTJets_V02.43_tree.root" ],
		[ "slimTTHadronic_V02.43_tree.root", "slimTTFull_V02.43_tree.root", "slimTTSemi_V02.43_tree.root"]
		)

	compareBinnedSamples(
		[ "slimZGammaNuNu3_400_inf_V03.00_tree.root" ],
		[ "slimZGamma_V02.43_tree.root" ],
		)

	inclusiveAndIsrSamples(
		[ "slimZGammaNuNu3_400_inf_V03.00_tree.root" ],
		[ "slimZGammaNuNu_V02.43_tree.root" ]
		)

	inclusiveAndIsrSamples(
		[ "slimZGammaNuNu3_400_inf_V03.00_tree.root" ],
		[ "slimZGammaNuNu_V02.43_tree.root" ],
		"pt130"
		)

	inclusiveAndIsrSamples(
		["slimTTJets_V02.43_tree.root"],
		["slimTTGamma_V02.43_tree.root" ]
		)

	inclusiveAndIsrSamples(
		[ "slimWJets_250_300_V02.43_tree.root", "slimWJets_300_400_V02.43_tree.root", "slimWJets_400_inf_V02.43_tree.root" ],
		["slimWGamma_50_130_V02.43_tree.root", "slimWGamma_130_inf_V02.43_tree.root" ]
		)

