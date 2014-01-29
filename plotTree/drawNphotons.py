#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from treeFunctions import *
Styles.tdrStyle2D()

ROOT.gStyle.SetPaintTextFormat("4.f")

def changeAxisTitle( h ):
	replacements = { "#gamma": "#ggamma", "#gamma_{jet}": "#fgamma", "#gamma_{e}": "#egamma" }

	for ax in h.GetYaxis(),h.GetXaxis(),h.GetZaxis():
		for old, new in replacements.iteritems():
			if ax.GetTitle() == old:
				ax.SetTitle( new )

def makeNphotonPlot( h3D, storeZero, datasetAbbr ):

	changeAxisTitle( h3D )


	paperWidth2 = 5.7
	ROOT.gStyle.SetPaperSize(paperWidth2,50.)
	infoText = ROOT.TLatex(0,.97, "#text{Private Work }#SI{19.8}{fb^{-1}}#, #sqrt{s}=#SI{8}{TeV}" )
	infoText.SetNDC()
	infoText.SetTextSize(1./9.81241)

	can = ROOT.TCanvas()
	can.cd()
	for combi in [ "yx", "yz", "zx" ]:
		h2 = h3D.Project3D( combi )
		if not storeZero:
			h2.SetBinContent( h2.FindBin(0,0), 0 )
		h2.SetMarkerSize(4)
		h2.GetYaxis().SetNdivisions(5,0,0)
		h2.GetXaxis().SetNdivisions(5,0,0)

		h2.SetLabelSize( 1./9.81241, "xyz" )
		h2.SetTitleSize( 1./9.81241, "xy" )
		h2.SetTitleOffset( 0.8, "x" )
		h2.SetTitleOffset( 0.7, "y" )
		h2.SetLabelOffset(0.0, "z" )


		h2.Draw("colz text")
		combi = combi.replace("x","g").replace("y","f").replace("z","e")
		infoText.Draw()
		SavePad( "nPhotons_%s_%s"%(datasetAbbr,combi) )
		SavePad( "nPhotons_%s_%s"%(datasetAbbr,combi), "~/master/documents/thesis/plots", ["tex"] )
		correctTiksPlot( "/home/knut/master/documents/thesis/plots/nPhotons_%s_%s.tex"%(datasetAbbr,combi))

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--zero", action="store_true" )
	arguments.add_argument("--add", action="store_true" )
	opts = arguments.parse_args()

	if opts.add:
		h3D = readHisto( opts.filenames[0], "nPhotons" )
		for filename in opts.filenames[1:]:
			h3D.Add( readHisto ( filename, "nPhotons" ) )
		makeNphotonPlot( h3D, opts.zero, shortName( opts.filenames) )
	else:
		for inName in opts.filenames:
			h3D = readHisto( inName, "nPhotons" )
			makeNphotonPlot( h3D, opts.zero, getDatasetAbbr(inName) )

