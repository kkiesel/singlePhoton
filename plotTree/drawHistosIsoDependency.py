#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from multiplot import Multihisto
from treeFunctions import *

def split2Din1DMultihist( h2D, useYaxis, uFlow, oFlow ):
	"""Creates many TH1 out of one TH2 and put them in a Multihisto.

	h2D: TH2 (rebin appropriate if necessary )
	useYaxis: The TH2 will be splitted in different values of Y.
		If false, the X axis will be used (bool).
	uFlow: Take the underFlow into account (bool).
	oFlow: Take the overFlow into account (bool).
	"""
	mh = Multihisto()
	style = ROOT.gROOT.GetStyle("tdrStyle") # for nice color gradient

	if useYaxis:
		axis = h2D.GetYaxis()
	else:
		axis = h2D.GetXaxis()

	variable = axis.GetTitle()
	nBins = axis.GetNbins()

	binRange = range( nBins + 2 )
	if not uFlow:
		binRange = binRange[1:]
	if not oFlow:
		binRange = binRange[:-1]

	for iColor, iBin in enumerate(binRange):
		if useYaxis:
			h = h2D.ProjectionX( "_px%s"%iBin, iBin, iBin )
		else:
			h = h2D.ProjectionY( "_px%s"%iBin, iBin, iBin )
		lowEdge = axis.GetBinLowEdge( iBin )
		upEdge = axis.GetBinUpEdge( iBin )

		if iBin == 0:
			title = "       %s < %s"%( variable, upEdge )
		elif iBin != nBins+1:
			title = "%s #leq %s < %s"%( lowEdge, variable, upEdge )
		else:
			title = "%s #leq %s"%( lowEdge, variable )

		color = style.GetColorPalette( int( 1.*iColor/(len(binRange)-1)*(style.GetNumberOfColors()-1) ) )
		h.SetLineColor( color )
		h.SetMarkerColor( h.GetLineColor() )
		h.SetMarkerSize(0)
		h.GetXaxis().SetTitle( h2D.GetXaxis().GetTitle() )
		mh.addHisto( h, title,draw="e0" )

	return mh

def drawDependency( inName, histName, rebinX, rebinY, axis="y", uFlow=True, oFlow=True ):
	hist2D = readHisto( inName, histName )
	hist2D = rebin2D( hist2D, rebinX, rebinY )

	mh = split2Din1DMultihist( hist2D, axis, uFlow, oFlow )
	for h in mh.histos:
		if h[0].Integral():
			h[0].Scale(1./h[0].Integral("width"),"width")
			print h[0].Integral("width")
		h[0].GetYaxis().SetTitle("Normed Entries, divided by bin width")

	mh.leg.SetX1(.5)
	if axis == "x" and "sigma" in hist2D.GetYaxis().GetTitle():
		mh.leg.SetY1(0.2)
		mh.leg.SetY2(0.6)

	can = ROOT.TCanvas("title", "name", 1000, 1400 )
	mh.Draw()
	datasetAbbr = getDatasetAbbr( inName )
	label = createDatasetLabel( datasetAbbr )
	label.Draw()
	SaveAs( can, "2Dto1D_%s_%s_%s_%s%s.pdf"%(histName,datasetAbbr,axis,len(rebinX),len(rebinY)) )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+",type=isValidFile )
	opts = arguments.parse_args()

	metBinning = readAxisConf( "met" )[2]
	metBinningCompact = [x for x in metBinning if x>=100 and x<500]
	photonBinning = [80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 600, 650,  700 ]

	for inName in opts.filenames:
		drawDependency( inName, "metSigma", metBinning, [0.0, 0.012, 0.014 ], True, False, False )
		drawDependency( inName, "metSigma", metBinning, [0.001*x for x in range(0,15,2) ], True, False )

		drawDependency( inName, "metChIso", metBinning, [ 0., 2.6, 15 ], True , False, False )
		drawDependency( inName, "metChIso", metBinning, [ 0, 0.5, 2,4,8,12,16, 22 ], True, False )

		drawDependency( inName, "metNeIso", metBinning, [ 3.5, 10 ], True, True )
		drawDependency( inName, "metNeIso", metBinning, range(0,30, 3), True, True )

		drawDependency( inName, "metPhIso", metBinning, [ 1.3, 10 ], True, True  )
		drawDependency( inName, "metPhIso", metBinning, range(0,30, 3), True, True )

		drawDependency( inName, "metHE", metBinning, [ 0., 0.05, 0.1 ], True, False  )
		drawDependency( inName, "metHE", metBinning, [0.05*x for x in range(13)], True, False )

