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


def getS7PuHisto():
	# Distribution used for S7 Summer2012 MC.

	Summer2012_S7 = [
		2.344E-05,
		2.344E-05,
		2.344E-05,
		2.344E-05,
		4.687E-04,
		4.687E-04,
		7.032E-04,
		9.414E-04,
		1.234E-03,
		1.603E-03,
		2.464E-03,
		3.250E-03,
		5.021E-03,
		6.644E-03,
		8.502E-03,
		1.121E-02,
		1.518E-02,
		2.033E-02,
		2.608E-02,
		3.171E-02,
		3.667E-02,
		4.060E-02,
		4.338E-02,
		4.520E-02,
		4.641E-02,
		4.735E-02,
		4.816E-02,
		4.881E-02,
		4.917E-02,
		4.909E-02,
		4.842E-02,
		4.707E-02,
		4.501E-02,
		4.228E-02,
		3.896E-02,
		3.521E-02,
		3.118E-02,
		2.702E-02,
		2.287E-02,
		1.885E-02,
		1.508E-02,
		1.166E-02,
		8.673E-03,
		6.190E-03,
		4.222E-03,
		2.746E-03,
		1.698E-03,
		9.971E-04,
		5.549E-04,
		2.924E-04,
		1.457E-04,
		6.864E-05,
		3.054E-05,
		1.282E-05,
		5.081E-06,
		1.898E-06,
		6.688E-07,
		2.221E-07,
		6.947E-08,
		2.047E-08 ]

	h = ROOT.TH1F(randomName(), "", 60, 0, 60 )
	for bin, content in enumerate( Summer2012_S7 ):
		h.SetBinContent( bin+1, content )

	return h


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	#arguments.add_argument("filenames", nargs="+", type=isValidFile )
	opts = arguments.parse_args()

	#hSig = getHistosFromFiles( opts.filenames )
	#hSig.SetName("pileup")
	#hSig.SetLineColor(2)

	hMc = readHisto( "../TreeWriter/pileUpReweighting/nTrueVertexMC.root", "pileupScenarioS10" )
	hSig = hMc.Clone()

	hData = readHisto( "../TreeWriter/pileUpReweighting/nTrueVertexData.root", "pileup" )
	hData.SetLineColor(1)
	hData.SetMarkerStyle(20)
	hData.Sumw2()

	hS7 = getS7PuHisto()
	hS7.SetLineColor(3)

	for h in hSig, hMc, hData, hS7:
		h.SetTitle(";true number vertices;[a.u.]")
		h.Scale( 1./h.Integral() )

	scaledTo = hS7.Clone( randomName() )
	scaledTo.Divide( hMc )
	scaledTo.Multiply( hData )
	scaledTo.SetLineColor(7)
	print scaledTo.Integral()


	mh = Multihisto()
	mh.setMaximum(0.07)
	mh.addHisto( hData, "Data", draw="p" )
	mh.addHisto( hMc, "S10 Scenario", draw="hist" )
	mh.addHisto( hSig, "Signal", draw="hist ")
	mh.addHisto( hS7, "S7 Scenario", draw="hist")
	mh.addHisto( scaledTo, "S7 Scenario (wrong scaled)", draw="hist")
	mh.Draw()

	ROOT.gPad.SetLogy(0)
	ROOT.gPad.SaveAs("plots/pileupScenarios.pdf")

