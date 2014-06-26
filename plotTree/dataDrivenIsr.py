#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from predictions import *
from qcdClosure import drawWeightHisto

ROOT.gStyle.SetOptLogy(0)

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

	ROOT.gPad.SaveAs("plots/chi2Minimization_%s.pdf"%plot)


def getkFactor( dataFiles, bkgFiles, plot, cut ):

	kFactorMin = 1
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

	bkgFiles = []
	bkgFiles.append( "slimZGammaLL_V02.%s_tree.root"%treeVersion )
	bkgFiles.append( "slimTTGamma_V03.%s_tree.root"%treeVersion )
	#bkgFiles.extend( ["slimGJets_400_inf_V03.%s_tree.root"%treeVersion, "slimGJets_200_400_V03.%s_tree.root"%treeVersion ] )

	kFactor = getkFactor( dataFiles, bkgFiles, opts.plot, chi2Cut )

	zgammall = getHists( ["slimZGammaLL_V02.%s_tree.root"%treeVersion], opts.plot, cut )
	zgammall.SetLineColor(2)
	zgammaIntegral,zgammaIntegralError = integralAndError(zgammall, zgammall.FindBin(60), zgammall.FindBin(119 ), "width" )

	ttgamma = getHists( ["slimTTGamma_V03.%s_tree.root"%treeVersion], opts.plot, cut )
	ttgamma.SetLineColor(4)

	#gjets = getHists( ["slimGJets_400_inf_V03.%s_tree.root"%treeVersion, "slimGJets_200_400_V03.%s_tree.root"%treeVersion ], opts.plot, cut )
	#gjets.SetLineColor( ROOT.kCyan )

	#qcd = getHists( ["slimQCD_1000_inf_V03.22_tree.root", "slimQCD_250_500_V03.22_tree.root", "slimQCD_500_1000_V03.22_tree.root"], opts.plot, cut )
	#qcd.SetLineColor( ROOT.kCyan+3 )

	#wjets = getHists( ["slimWJets_250_300_V03.22_tree.root", "slimWJets_300_400_V03.22_tree.root", "slimWJets_400_inf_V03.22_tree.root" ], opts.plot, cut )
	#wjets.SetLineColor( ROOT.kGreen+4 )

	#wgamma = getHists( ["slimWGamma_130_inf_V03.22_tree.root", "slimWGamma_50_130_V03.22_tree.root" ], opts.plot, cut )
	#wgamma.SetLineColor( ROOT.kGreen-4 )

	#signal = getHists( ["slimW_1700_720_375_V03.22_tree.root" ], opts.plot, cut )
	#signal.SetLineColor( ROOT.kGreen )
	#signal.SetLineWidth(2)


	for h in zgammall,ttgamma:#, gjets:
		h.Scale( kFactor )

	mh = Multihisto()
	mh.setMinimum(0)
	mh.addHisto( data, "Data", draw="pe" )
	mh.addHisto( zgammall, "#gammaZ(ll)", True )
	mh.addHisto( ttgamma, "#gammat#bar{t}", True )
	#mh.addHisto( gjets, "#gammaJet", True )
	#mh.addHisto( qcd, "Multijet", True )
	#mh.addHisto( wjets, "W", True )
	#mh.addHisto( wgamma, "#gammaW", True )
	#mh.addHisto( signal, "Wino", False )

	mh.Draw()

	fitFunc = ROOT.TF1( "fitfunc", "gaus(0)+pol1(3)", 20, 140 )
	fitFunc.SetParameters( 2, 91, 13, -0.1, 0.8 )
	fitFunc.FixParameter(1, 91 )
	data.Fit( "fitfunc", "IMRN" )
	fitFunc.SetLineColor(1)
	fitFunc.SetLineWidth(2)
	fitFunc.Draw("same")
	fitDataIntegral =  fitFunc.Integral(60, 120)
	fitDataIntegral -= ( fitFunc.GetParameter(3)*(120-60) + fitFunc.GetParameter(4)*(120**2-60**2)/2 ) # subtract integral of pol1
	fitDataIntegralError =  fitFunc.IntegralError(60, 120)

	#gausFunc = ROOT.TF1( "gausfunc", "gaus", 60, 120 )
	#gausFunc.SetParameters( fitFunc.GetParameter(0), fitFunc.GetParameter(1), fitFunc.GetParameter(2) )
	#gausFunc.SetParError( 0, fitFunc.GetParError(0) )
	#gausFunc.SetParError( 1, fitFunc.GetParError(1) )
	#gausFunc.SetParError( 2, fitFunc.GetParError(2) )
	#fitDataIntegral =  gausFunc.Integral(60, 120)
	#fitDataIntegralError =  gausFunc.IntegralError(60, 120)
	kFactorFit = fitDataIntegral/zgammaIntegral
	kFactorFitError = kFactorFit * sqrt( (zgammaIntegralError / zgammaIntegral)**2 + (fitDataIntegralError/fitDataIntegral)**2 )

	print "k-factor(Fit) = %s ± %s"%(kFactorFit, kFactorFitError )


	ROOT.gPad.SaveAs( "plots/isrkFactor_%s.pdf"%opts.plot )

