#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from treeFunctions import *
Styles.tdrStyle2D()

def makeNphotonPlot( h3D, storeZero, datasetAbbr ):

	can = ROOT.TCanvas()
	can.cd()
	for combi in [ "yx", "yz", "zx" ]:
		h2 = h3D.Project3D( combi )
		h2.GetXaxis().SetNdivisions(8,0,0)
		h2.GetYaxis().SetNdivisions(8,0,0)
		if not storeZero:
			h2.SetBinContent( h2.FindBin(0,0), 0 )
		h2.Draw("colz text")
		can.SaveAs( "plots/nPhotons_%s_%s.pdf"%(datasetAbbr,combi) )

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
		makeNphotonPlot( h3D, opts.zero, "sum" )
	else:
		for inName in opts.filenames:
			h3D = readHisto( inName, "nPhotons" )
			makeNphotonPlot( h3D, opts.zero, getDatasetAbbr(inName) )

