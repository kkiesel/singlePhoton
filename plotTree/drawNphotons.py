#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from treeFunctions import *

def makeNphotonPlot( inName ):
	datasetAbbr = getDatasetAbbr( inName )
	h3D = readHisto( inName, "nPhotons" )

	can = ROOT.TCanvas()
	can.cd()
	for combi in [ "yx", "yz", "xz" ]:
		h2 = h3D.Project3D( combi )
		h2.GetXaxis().SetNdivisions(8,0,0)
		h2.GetYaxis().SetNdivisions(8,0,0)
		h2.Draw("colz text")
		can.SaveAs( "plots/nPhotons_%s_%s.pdf"%(datasetAbbr,combi) )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	opts = arguments.parse_args()

	for inName in opts.filenames:
		makeNphotonPlot( inName )

