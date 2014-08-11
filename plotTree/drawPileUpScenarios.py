#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from predictions import *


def getHistosFromFiles( filenames, histname="trueNVtx" ):
	h = readHisto( filenames[0], histname )
	for filename in filenames[1:]:
		h.Add( readHisto( filename, histname ) )

	return h

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	#arguments.add_argument("--filenames", nargs="+", type=isValidFile )
	opts = arguments.parse_args()

	s10 = readHisto( "../TreeWriter/pileUpReweighting/nTrueVertexMC.root", "pileupScenarioS10" )
	s7 = readHisto( "../TreeWriter/pileUpReweighting/nTrueVertexMC.root", "pileupScenarioS7" )
	data = readHisto( "../TreeWriter/pileUpReweighting/nTrueVertexData.root", "pileup" )

	data.SetLineColor(1)
	data.SetMarkerStyle(20)

	s10.SetLineColor(2)
	s7.SetLineColor( ROOT.kBlue )

	for h in data, s7, s10:
		h.SetTitle(";true number of vertices;[a.u.]")
		h.Scale( 1./h.Integral() )
		h.GetXaxis().SetRangeUser( 5, 50 )

	S10scaledUsingS10 = s10.Clone( randomName() )
	S10scaledUsingS10.Divide( s10 )
	S10scaledUsingS10.Multiply( data )
	S10scaledUsingS10.SetLineStyle( 2 )

	S7scaledUsingS7 = s7.Clone( randomName() )
	S7scaledUsingS7.Divide( s7 )
	S7scaledUsingS7.Multiply( data )
	S7scaledUsingS7.SetLineStyle( 3 )

	S7scaledUsingS10 = s7.Clone( randomName() )
	S7scaledUsingS10.Divide( s10 )
	S7scaledUsingS10.Multiply( data )
	S7scaledUsingS10.SetLineStyle( 8 )

	mh = Multihisto()
	mh.setMaximum(0.07)
	mh.addHisto( data, "Data", draw="p" )
	mh.addHisto( s10, "S10 Scenario", draw="hist" )
	mh.addHisto( S10scaledUsingS10, "S10 scaled by data/S10", draw="hist")
	mh.addHisto( s7, "S7 Scenario", draw="hist")
	mh.addHisto( S7scaledUsingS7, "S7 scaled by data/S7", draw="hist")
	mh.addHisto( S7scaledUsingS10, "S7 scaled by data/S10", draw="hist")
	mh.Draw()
	mh.leg.SetX1( .65 )
	mh.leg.SetY1( .7 )
	mh.leg.SetX2( 1 )
	mh.leg.SetY2( 1 )

	ROOT.gPad.SetLogy(0)
	ROOT.gPad.SaveAs("plots/pileupScenarios.pdf")

	for name, h in [("data",data), ("s7", s7), ("s10", s10), ("s10scaleds10", S10scaledUsingS10), ("s7sceleds7", S7scaledUsingS7), ("s7sceleds10", S7scaledUsingS10)]:
		print name, h.GetBinContent(h.FindBin(20))


##################### New plot ###################################
	# pu output distributions
	signalS7 = readHisto( "puInputs/testS7.root", "nTrueVertex" )
	signalS10 = readHisto( "puInputs/testS10.root", "nTrueVertex" )
	signalS0 = readHisto( "puInputs/testNoWeight.root", "nTrueVertex" )

	for name, hist in [ ("S7 ", signalS7), ("S10", signalS10), ("NPu", signalS0 ) ]:
		print name, hist.Integral(0, hist.FindBin(20)-1 )/600, hist.Integral(hist.FindBin(20),hist.FindBin(30)-1)/600, hist.Integral(hist.FindBin(30), -1 )/600

	signalS7.SetLineColor(1)
	signalS10.SetLineColor(2)

	for h in signalS7, signalS10, signalS0:
		h.SetTitle(";true number of vertices;Events")
		h.GetXaxis().SetRangeUser( 5, 50 )


	mh2 = Multihisto()
	mh2.leg.SetX1( .65 )
	mh2.leg.SetY1( .7 )
	mh2.leg.SetX2( 1 )
	mh2.leg.SetY2( 1 )

	mh2.addHisto( signalS7, "weighted by data/S7", draw="e" )
	mh2.addHisto( signalS10, "weighted by data/S10", draw="e" )
	mh2.addHisto( signalS0, "no pu weighting", draw="e" )
	mh2.Draw()

	ROOT.gPad.SaveAs("plots/pileupDifferentWeightingSignal.pdf")

