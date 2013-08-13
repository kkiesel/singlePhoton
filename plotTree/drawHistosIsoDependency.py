#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import ROOT
import argparse
import ConfigParser
import re
from math import sqrt
from multiplot import *
from treeFunctions import *
import ratios
import Styles

Styles.tdrStyle()
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libTreeObjects.so")

# to use user defined help message, sys.arv has to be sent to python and not
# to TApplication
ROOT.PyConfig.IgnoreCommandLineOptions = True

def split2Din1DMultihist( h2D, axis, uFlow, oFlow ):
	mh = Multihisto()
	style = ROOT.gROOT.GetStyle("tdrStyle") # for nice color gradient

	if axis == "y":
		variable = h2D.GetYaxis().GetTitle()
		binRange = range( h2D.GetNbinsY()+2 )
		if not uFlow:
			binRange = binRange[1:]
		if not oFlow:
			binRange = binRange[:-1]
		for yBin in binRange:
			h = h2D.ProjectionX( "_px%s"%yBin, yBin, yBin )
			if yBin == 0:
				title = "       %s < %s"%(variable, h2D.GetYaxis().GetBinUpEdge(yBin))
			elif yBin != h2D.GetNbinsY()+1:
				title = "%s #leq %s < %s"%(h2D.GetYaxis().GetBinLowEdge(yBin), variable, h2D.GetYaxis().GetBinUpEdge(yBin) )
			else:
				title = "%s #leq %s"%(h2D.GetYaxis().GetBinLowEdge(yBin),variable )

			h.SetLineColor( style.GetColorPalette(int(1.*yBin/(len(binRange))*(style.GetNumberOfColors()-1) ) ) )
			h.SetMarkerColor( h.GetLineColor() )
			h.SetMarkerSize(0)
			h.GetXaxis().SetTitle( h2D.GetXaxis().GetTitle() )
			mh.addHisto( h, title,draw="e0" )
	if axis == "x":
		variable = h2D.GetXaxis().GetTitle()
		binRange = range( h2D.GetNbinsX()+2 )
		if not uFlow:
			binRange = binRange[1:]
		if not oFlow:
			binRange = binRange[:-1]
		for yBin in binRange:
			h = h2D.ProjectionY( "_px%s"%yBin, yBin, yBin )
			if yBin == 0:
				title = "       %s < %s"%(variable, h2D.GetXaxis().GetBinUpEdge(yBin))
			elif yBin != h2D.GetNbinsX()+1:
				title = "%s #leq %s < %s"%(h2D.GetXaxis().GetBinLowEdge(yBin), variable, h2D.GetXaxis().GetBinUpEdge(yBin) )
			else:
				title = "%s #leq %s"%(h2D.GetXaxis().GetBinLowEdge(yBin),variable )

			h.SetLineColor( style.GetColorPalette(int(1.*yBin/(len(binRange)-1)*(style.GetNumberOfColors()-1) ) ) )
			h.SetMarkerColor( h.GetLineColor() )
			h.SetMarkerSize(0)
			h.GetXaxis().SetTitle( h2D.GetYaxis().GetTitle() )
			mh.addHisto( h, title,draw="e0" )

	return mh

def drawDependency( inName, histName, rebinX, rebinY, axis="y", uFlow=True, oFlow=True ):
	hist2D = readHisto( inName, histName )
	hist2D = rebin2D( hist2D, rebinX, rebinY )

	mh = split2Din1DMultihist( hist2D, axis, uFlow, oFlow )
	for h in mh.histos:
		if h[0].Integral():
			h[0].Scale(1./h[0].Integral(),"width")
		h[0].GetYaxis().SetTitle("Normed Entries, divided by bin width")

	if axis == "x" and "sigma" in hist2D.GetYaxis().GetTitle():
		mh.leg.SetY1(0.2)
		mh.leg.SetY2(0.6)

	can = ROOT.TCanvas()
	mh.Draw()
	datasetAffix = re.match("(slim)?(.*)_V.*", inName ).groups()[-1]
	can.SaveAs("plots/2Dto1D_%s_%s_%s_%s%s.pdf"%(histName,datasetAffix,axis,len(rebinX),len(rebinY)))

if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Calculate weighting "
			+"factors for QCD background estimation." )
	arguments.add_argument("--input", default=[""], nargs="+" )
	opts = arguments.parse_args()

	metBinning = readAxisConf( "met" )[2]
	metBinningCompact = [x for x in metBinning if x>=100 and x<500]
	photonBinning = [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600, 650,  700 ]

	for inName in opts.input:
		drawDependency( inName, "metSigma", metBinning, [0.0, 0.012, 0.014 ], "y", False, False )
		drawDependency( inName, "metSigma", metBinning, [0.001*x for x in range(0,15,2) ] )

		drawDependency( inName, "metChIso", metBinning, [ 0., 2.6 ], "y", False  )
		drawDependency( inName, "metChIso", metBinning, [ 0, 0.5, 2,4,8,12,16, 22 ], "y", False )

		drawDependency( inName, "metNeIso", metBinning, [ 0., 2.6 ], "y", False  )
		drawDependency( inName, "metNeIso", metBinning, range(0,30, 3), "y", False )

		drawDependency( inName, "metPhIso", metBinning, [ 0., 2.6 ], "y", False  )
		drawDependency( inName, "metPhIso", metBinning, range(0,30, 3), "y", False )


		#drawDependency( inName, "metSigma", metBinningCompact, [ 0.001*x for x in range(16)], "x" )
		#drawDependency( inName, "ptSigma", [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600 ], [ 0.006, 0.012, 0.014 ] )
		#drawDependency( inName, "ptJetSigma", [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600 ], [ 0.006, 0.012, 0.014 ] )
		#drawDependency( inName, "ptSigma", [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600 ], [ 0.001*x for x in range(16)], "x" )
		#drawDependency( inName, "ptJetSigma", [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600 ], [ 0.001*x for x in range(16)], "x" )
		#drawDependency( inName, "ptChIso", [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600, 1000, 1500 ], [ 0.7, 2.6, 4.5, 10 ] )
		#drawDependency( inName, "ptJetChIso", [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600, 1000,1500 ], [ 0.7, 2.6, 4.5, 10 ] )
		#drawDependency( inName, "ptJetChIso", [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600 ], [ 0.7, 2.6, 4.5, 10 ], "x" )
		#drawDependency( inName, "htSigma", range( 400, 1100, 50), [ 0.006, 0.012, 0.014 ] )
		#drawDependency( inName, "htSigma", [300, 500, 600, 700, 800,900 ], [ 0.001*x for x in range(16)], "x" )
		#drawDependency( inName, "metPt", metBinning, photonBinning )

