#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

def printBin( histo, error=False, scale=1. ):
	out = ""
	hClone = histo.Clone( randomName() )
	hClone.Scale( scale )
	for bin in range(1, hClone.GetNbinsX()+2):
		val = hClone.GetBinError(bin) if error else hClone.GetBinContent(bin)
		out += " %s"%val
	return out


def createSignalCard( filename, binningStr ):

	# binning
	binList = readAxisConf( "met" )[2]
	binList = [x for x in binList if x>=100] # only values above 100
	binList = binList[0:-1] # exclude last bin

	binListDict = {}
	binListDict['optimized'] = [ 100, 270, 350, 500 ]
	binListDict['old'] = [ 100, 120, 160, 200, 270, 350, 500 ]
	binListDict['new'] = [ 100, 120, 150, 200, 260, 340, 430, 550 ]
	binListDict['fibo'] = [ 100, 170, 285, 455, 740 ]

	binList = binListDict[binningStr]
	import array
	metBinning = array.array( "d", binList )

	# christian wants to have a global counter
	# no physical information
	counter = 0

	fakeRateSysError = 11./100

	import time
	signalCardString = """
# This file contains event yield and uncertainties for signal points.
# Produced by:
#    Knut Kiesel
# on:
#    %s
# with input file:
#    %s
# There will be a common section for all signal points, folllowed by the information
# for each signal point. The numbering of the signal points is arbitrary.
# For bin edges, the ROOT conversion is used (lower bin edge included, upper bin
# edge excluded).
"""%( time.strftime("%Y-%m-%d %H:%M:%S"), filename )

	scanName = "unknown"
	mBino = 0
	mWino = 0
	nGen = 0
	if "W_gsq_" in filename:
		scanName = "Spectra_gsq_W"
		mWino = 375
		mBino = 5000
		nGen = 60000
	if "B_gsq_" in filename:
		scanName = "Spectra_gsq_B"
		mBino = 375
		mWino = 5000
		nGen = 10000

	commonInformation = """
signal scan = %s
bino mass = %s
wino mass = %s
nGen = %s
nMetBins = %s
""" %(scanName, mBino, mWino, nGen, len(binList))

	binListInf = binList + ['infinity']

	for i in range(len(binList)):
		commonInformation += "bin %s= %s to %s\n"%(i,binListInf[i],binListInf[i+1])

	separationLine = "\n###############################################################\n"

	signalCardString += separationLine
	signalCardString += commonInformation
	signalCardString += separationLine

	# loop over signal scans
	xList = range( 400, 2001, 100 )
	yList = range( 420, 2021, 100 )

	for x in xList:
		for y in yList:

			signalCardString += separationLine
			if "W_gsq_" in filename:
				if x == 500 and y == 520 or x == 1000 and y == 1520:
					signalCardString += "Point not avaiable\n"
					continue
			if "B_gsq_" in filename:
				if x == 1900 and y == 1820:
					signalCardString += "Point not avaiable\n"
					continue

			signalCardString += "Point %s squark mass = %s\n"%(counter, x)
			signalCardString += "Point %s gluino mass = %s\n"%(counter, y)

			gHisto = readHisto( filename, "gMet%s_%s"%(x,y) ).Rebin( len(metBinning)-1, randomName(), metBinning )
			eHisto = readHisto( filename, "eMet%s_%s"%(x,y) ).Rebin( len(metBinning)-1, randomName(), metBinning )
			fHisto = readHisto( filename, "fMet%s_%s"%(x,y) ).Rebin( len(metBinning)-1, randomName(), metBinning )
			fMetError = readHisto( filename, "fMetError%s_%s"%(x,y) ).Rebin( len(metBinning)-1, randomName(), metBinning )
			gMetPuUp = readHisto( filename, "gMetPuUp%s_%s"%(x,y) ).Rebin( len(metBinning)-1, randomName(), metBinning )
			gMetPuDown = readHisto( filename, "gMetPuDown%s_%s"%(x,y) ).Rebin( len(metBinning)-1, randomName(), metBinning )
			gMetPu = gMetPuUp - gMetPuDown
			gMetPu.Scale(0.5)

			# assume the same binning for all histograms
			signalCardString += "Point %s number of signal events in bins = %s\n"%(counter, printBin( gHisto ))
			signalCardString += "Point %s statistical error of signal events in bins = %s\n"%(counter, printBin( gHisto, True ))
			signalCardString += "Point %s pu uncert = %s\n"%(counter, printBin( gMetPu ))
			signalCardString += "Point %s EWK prediction = %s\n"%(counter, printBin( eHisto ))
			signalCardString += "Point %s EWK statistical error from signal = %s\n"%(counter, printBin( eHisto, True ))
			signalCardString += "Point %s EWK statistical error from yutaro = %s\n"%(counter, printBin( eHisto, scale=fakeRateSysError ))
			signalCardString += "Point %s QCD prediction = %s\n"%(counter, printBin( fHisto ) )
			signalCardString += "Point %s QCD statistical error from signal = %s\n"%(counter, printBin( fHisto, True ))
			signalCardString += "Point %s QCD systematical error from weight = %s\n"%(counter, printBin( fMetError ) )

			counter += 1

	#scanName += "" if len(binList) == 6 else "_%smetBins"%len(binList)
	#scanName += binningStr

	outputFileName = "eventYield%s-%s.txt"%(scanName, time.strftime("%Y-%m-%d"))
	signalCardFile = open( outputFileName, "w")
	signalCardFile.write( signalCardString )
	signalCardFile.close()
	print "Wrote to %s"%outputFileName


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile, default= "W_gsq_V02.45.root" )
	opts = arguments.parse_args()

	for filename in opts.filenames:
		#for binningStr in 'optimized', 'old', 'new', 'fibo':
		for binningStr in ['old']:
			createSignalCard( filename, binningStr )




