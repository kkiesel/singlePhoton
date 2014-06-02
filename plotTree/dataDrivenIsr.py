#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from predictions import *
from qcdClosure import drawWeightHisto

ROOT.gStyle.SetOptLogy(0)

def drawChi2( tuples ):

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

	ROOT.gPad.SaveAs("plots/chi2MinimizationZmumu.pdf")
	ROOT.gPad.SaveAs("plots/chi2MinimizationZmumu.png")


def getkFactor( dataFiles, bkgFiles, cut ):

	kFactorMin = 1
	kFactorMax = 4
	kFactorN = 100

	data = getHists( dataFiles, "mll", cut )
	bkg = getHists( bkgFiles, "mll", cut )

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
		chi2 /= ndf

		#chi2 = data.Chi2Test( scaledbkg, "uw CHI2/NDF p" )
		#ROOT.gPad.SaveAs("test%s.pdf"%i)
		#print kFactor, chi2
		tuples.append( (kFactor, chi2) )

	drawChi2( tuples )

	return min( tuples, key=lambda t: t[1] )[0]


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("--plot", nargs="+", default = ["met"] )
	opts = arguments.parse_args()

	cut = "std::abs(photons[0].eta)<1.4442 && mll>10 && met<100"
	chi2Cut = cut + "&& 60<mll && mll<120"

	treeVersion = 24
	dataFiles = [ "PhotonHad%s_V03.%s_tree.root"%(x,treeVersion) for x in ["A","B","C", "D" ] ]
	data = getHists( dataFiles, "mll", cut )

	bkgFiles = []
	bkgFiles.append( "slimZGammaLL_V02.24_tree.root" )
	bkgFiles.append( "slimTTGamma_V03.24_tree.root" )
	#bkgFiles.extend( ["slimGJets_400_inf_V03.22_tree.root", "slimGJets_200_400_V03.22_tree.root" ] )

	kFactor = getkFactor( dataFiles, bkgFiles, chi2Cut)

	h2 = getHists( ["slimZGammaLL_V02.24_tree.root"], "mll", cut )
	h2.SetLineColor(2)

	h3 = getHists( ["slimTTGamma_V03.24_tree.root"], "mll", cut )
	h3.SetLineColor(4)

	#gjets = getHists( ["slimGJets_400_inf_V03.22_tree.root", "slimGJets_200_400_V03.22_tree.root" ], "mll", cut )
	#gjets.SetLineColor( ROOT.kCyan )

	#qcd = getHists( ["slimQCD_1000_inf_V03.22_tree.root", "slimQCD_250_500_V03.22_tree.root", "slimQCD_500_1000_V03.22_tree.root"], "mll", cut )
	#qcd.SetLineColor( ROOT.kCyan+3 )

	#wjets = getHists( ["slimWJets_250_300_V03.22_tree.root", "slimWJets_300_400_V03.22_tree.root", "slimWJets_400_inf_V03.22_tree.root" ], "mll", cut )
	#wjets.SetLineColor( ROOT.kGreen+4 )

	#wgamma = getHists( ["slimWGamma_130_inf_V03.22_tree.root", "slimWGamma_50_130_V03.22_tree.root" ], "mll", cut )
	#wgamma.SetLineColor( ROOT.kGreen-4 )

	signal = getHists( ["slimW_1700_720_375_V03.22_tree.root" ], "mll", cut )
	signal.SetLineColor( ROOT.kGreen )
	signal.SetLineWidth(2)


	for h in h2,h3:
		h.Scale( kFactor )

	mh = Multihisto()
	mh.setMinimum(0)
	mh.addHisto( data, "Data", draw="pe" )
	mh.addHisto( h2, "#gammaZ(ll)", True )
	mh.addHisto( h3, "#gammat#bar{t}", True )
	#mh.addHisto( gjets, "#gammaJet", True )
	#mh.addHisto( qcd, "Multijet", True )
	#mh.addHisto( wjets, "W", True )
	#mh.addHisto( wgamma, "#gammaW", True )
	#mh.addHisto( signal, "Wino", False )

	mh.Draw()

	ROOT.gPad.SaveAs( "plots/data_invariantMassMuMu.pdf" )
	ROOT.gPad.SaveAs( "plots/data_invariantMassMuMu.png" )

