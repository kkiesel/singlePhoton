#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *
ROOT.gStyle.SetErrorX(0)
ROOT.gStyle.SetEndErrorSize(0)

def doTriggerPlot( ht=True ):
	#filenamePart = "/home/knut/master/triggerPlots/4knut/mergedHistos_%sTrigger.root"
	filenamePart = "/user/kiesel/root-files/mergedHistos_%s.root"
	filenamePart2 = "HT" if ht else "Photon"
	plotName = "photonTrigPt" if ht else "hthlt"

	f = ROOT.TFile( filenamePart%filenamePart2 )
	nom = f.Get( "select_All_%s_Nominator/PreselCut_%s"%(filenamePart2,plotName) )
	den = f.Get( "select_All_%s_Denominator/PreselCut_%s"%(filenamePart2,plotName) )

	if not den or not nom:
		print "histogram not found"
		return

	eff = ROOT.TGraphAsymmErrors( nom, den )
	eff.SetMinimum(0)
	eff.SetMaximum(1.05)
	xTitle = "H_{T}" if not ht else "p_{T^{*}}"
	eff.SetTitle(";%s [GeV];Efficiency"%xTitle )
	eff.SetLineColor(1)
	eff.SetMarkerColor(1)
	eff.SetLineWidth(1)
	if not ht:
		#eff.GetXaxis().SetRangeUser( 300, 1000 )
		eff.GetXaxis().SetLimits( 300, 1000 )
	else:
		#eff.GetXaxis().SetRangeUser( 50, 500 )
		eff.GetXaxis().SetLimits( 50, 500 )



	can = ROOT.TCanvas( randomName(), "", 600, 600 )
	can.SetLogy( 0 )
	can.SetTopMargin(0.08)

	eff.Draw("aep")

	infoText = ROOT.TLatex()
	infoText.SetNDC()
	infoText.SetTextFont( eff.GetXaxis().GetLabelFont() )
	infoText.SetTextSize( eff.GetXaxis().GetLabelSize() )
	infoText.SetText( .01, .95, "CMS Preliminary          19.7fb^{-1} (8 TeV)" )
	infoText.Draw()

	# lines
	cutValue = 110 if ht else 500
	l = ROOT.TLine( cutValue, 0, cutValue, eff.GetMaximum() )
	l.SetLineStyle(2)
	l.SetLineColor(2)
	l.SetLineWidth(2)
	l.Draw()

	oneLine = ROOT.TLine( eff.GetXaxis().GetXmin(), 1, eff.GetXaxis().GetXmax(), 1 )
	oneLine.SetLineStyle(3)
	oneLine.SetLineWidth(2)
	oneLine.Draw()

	outName = "plots/triggerEfficiency" + filenamePart2
	ROOT.gPad.SaveAs( outName + ".pdf" )

doTriggerPlot()
doTriggerPlot( False )
