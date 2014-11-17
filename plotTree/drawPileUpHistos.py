#! /usr/bin/env python2
# -*- coding: utf-8 -*-

import ROOT
from multiplot import Multihisto
from treeFunctions import *

def texStyle():
	style = ROOT.TStyle("texStyle", "Style to generate PGF/TikZ")

	# Canvas
	style.SetCanvasColor( ROOT.kWhite )
	style.SetCanvasBorderMode(0)
	style.SetCanvasDefH(1000)
	style.SetCanvasDefW(2000)

	# Pad:
	style.SetPadBorderMode(0)

	# Margins:
	style.SetPadTopMargin(0.05)
	style.SetPadBottomMargin(0.13)
	style.SetPadLeftMargin(0.16)
	style.SetPadRightMargin(0.02)

	# axis titles
	style.SetTitleFont(42, "xyz")
	style.SetTitleSize(0.09, "xyz")
	style.SetLabelSize(0.09, "xyz")

	# Statistic box
	style.SetOptStat(0)

	# Other
	style.SetPalette(1)
	# x-value is approximation of text width, y value large
	style.SetPaperSize(14.65,50.)

	style.SetPadTickX(1)
	style.SetPadTickY(1)

	style.cd()
	return style

def drawPileUpHistos( saveTex=False ):

	inputHistPath = "../TreeWriter/pileUpReweighting/"

	data = readHisto( inputHistPath + "nTrueVertexData.root", "pileup" )
	data.SetLineColor(1)
	data.SetMarkerColor(1)
	data.SetMarkerStyle(20)
	data.Sumw2()
	if saveTex:
		data.SetMarkerSize(.8)
	else:
		data.SetMarkerSize(1.2)

	mc = readHisto( inputHistPath + "nTrueVertexMC.root", "pileupS10" )
	mc.SetLineColor(2)
	mc.SetMarkerColor(2)

	for h in [mc, data ]:
		h.Scale( 1./h.Integral() )
		h.SetTitle( ";Number of Pile-up Events;Normalized Entries" )
		if saveTex:
			h.SetLineWidth(2)
			h.SetTitleSize( 0.06311227345609463, "yx" )
			h.SetLabelSize( 0.06311227345609463, "yx" )

	muhist = Multihisto()
	muhist.addHisto( data, "Data", draw="ep" )
	muhist.addHisto( mc, "Simulation", draw="hist" )
	muhist.setMinimum(0)

	if saveTex:
		muhist.leg.SetTextSize(0.063112267888)
		texStyle()


	pc1 = ROOT.TLatex(0,.96, "CMS Private Work")
	pc2 = ROOT.TLatex( .76,.96, "\SI{19.8}{fb^{-1}} #sqrt{s}=\SI{8}{TeV}")

	masterPath = "~/master/documents/thesis/plots/"
	can = ROOT.TCanvas("","",2000,1000)
	can.cd()
	can.SetLogy(0)
	muhist.Draw()
	for pc in [pc1, pc2]:
		pc.SetNDC()
		if saveTex:
			pc.SetTextSize(0.06311227345609463)
		pc.Draw()


	ending = 'tex' if saveTex else 'pdf'

	ROOT.gPad.SaveAs(masterPath+"puDistribution.%s"%ending)

drawPileUpHistos( True )
#drawPileUpHistos( False )
