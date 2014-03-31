#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
metFilters = [
		"Passing all filters",
		"CSCBeamHalo",
		"HcalNoise",
		"EcalDeadCellTP",
		"EcalDeadCellBE",
		"TrackingFailure",
		"EEBadSC",
		"HcalLaserOccupancy",
		"HcalLaserEventList",
		"HcalLaserRECOUserStep",
		"EcalLaserCorr",
		"ManyStripClus53X",
		"TooManyStripClus53X",
		"LogErrorTooManyClusters",
		"LogErrorTooManyTripletsPairs",
		"LogErrorTooManySeeds",
		"EERingOfFire",
		"InconsistentMuon",
		"GreedyMuon"
	]

def printMetFilters( h, sort ):
	dataToPrint = []
	for i in range(h.GetNbinsX()):
		dataToPrint.append( [ metFilters[i], str(int(h.GetBinContent(i))), "%.3f"%(h.GetBinContent(i)/h.GetBinContent(0)*100) ] )

	if sort:
		dataToPrint = sorted( dataToPrint, key=lambda entry: int(entry[1]), reverse=True )

	dataToPrint.insert( 0, ["Filter", "#Events", "ratio [%]"] )
	col_width = max(len(word) for row in dataToPrint for word in row) + 2  # padding
	for row in dataToPrint:
		print "".join(word.rjust(col_width) for word in row)


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Print out filters" )
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--sum", action="store_true" )
	arguments.add_argument( "--sort", action="store_true" )
	opts = arguments.parse_args()

	if opts.sum:
		print " ".join(mergeDatasetAbbr([getDatasetAbbr(x) for x in opts.filenames ]))
		h = readHisto( opts.filenames[0], "metFilters" )
		for filename in opts.filenames[1:]:
			h.Add( readHisto( filename, "metFilters" ) )
		printMetFilters( h, opts.sort )

	else:
		for filename in opts.filenames:
			print getDatasetAbbr(filename)
			h = readHisto( filename, "metFilters" )
			printMetFilters( h, opts.sort )

