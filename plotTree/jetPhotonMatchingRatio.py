#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from sys import stdout
from treeFunctions import *
import re




def ratioFrom2D( h2 ):
	h1 = ROOT.TH1F( randomName(), "ratio from 2D histogram", 1000, 0, 20 )
	for xbin in range( h2.GetNbinsX()+1 ):
		for ybin in range( h2.GetNbinsY()+1 ):
			h1.Fill( h2.GetYaxis().GetBinCenter(ybin) / h2.GetXaxis().GetBinCenter(xbin), h2.GetBinContent( xbin, ybin ) )


	can = ROOT.TCanvas()
	can.cd()
	h1.Draw()
	SaveAs( can, "test" )


ratioFrom2D( readHisto( "slimTTJets_V02.32_tree.root", "matchJetPt" ) )
