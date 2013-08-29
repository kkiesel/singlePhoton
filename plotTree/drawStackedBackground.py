#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

def drawStackedBackground( plot, treeName, listOfFiles ):
	mh = Multihisto()
	mh.orderByIntegral = False
	for iColor, fileName in enumerate(listOfFiles):
		datasetAbbr = getDatasetAbbr( fileName )
		tree = readTree( fileName, treeName )
		histo = getHisto( tree, plot, color=iColor+3 )
		mh.addHisto( histo, datasetAbbr, toStack=True, draw="hist" )

	can = ROOT.TCanvas()
	mh.Draw()
	can.SaveAs("plots/stackedHisto_%s_%s.pdf"%(treeName, plot))

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--tree", default="photonTree" )
	opts = arguments.parse_args()

	drawStackedBackground( opts.plot, opts.tree, opts.filenames )
