#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

def getFromEvent( event ):
	ht = event.photons[0].ptJet()
	for jet in event.jets:
		if jet.isMatch(3):
			ht += jet.pt
	return ht


def compareTrees( plot="photons.pt", filename="slimAllQCD_V02.28_tree.root" ):
	label, unit, binning = readAxisConf( "ht" )
	import array
	fH = ROOT.TH1F( randomName(), ";%s;"%label, len(binning)-1, array.array('d', binning) )
	fH.SetLineColor(2)
	fH.SetMarkerColor(2)
	fH.Sumw2()
	#cut = "!photons[0].isGen(0)"
	cut = "1"

	gTree = readTree( filename, "photonTree" )
	fTree = readTree( filename, "photonJetTree" )

	gH = getHisto( gTree, "ht", color=1, cut=cut,firstBin=0,lastBin=1400 )
	#fH = getHisto( fTree, plot, color=2, cut=cut,firstBin=0,lastBin=1400 )

	for event in fTree:
		fH.Fill( getFromEvent( event ) )


	for h in [gH, fH]:
		h.Scale( 1./h.Integral() )

	mh = Multihisto()
	mh.addHisto( gH, "new", draw="hist e" )
	mh.addHisto( fH, "old (recalculated)", draw="hist e" )

	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()


	from myRatio import Ratio
	r = Ratio( "#gamma/#gamma_{jet}", gH, fH )
	r.draw(0,2)


	SaveAs( can, "compare_%s_norm"%plot )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	#arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--tree", default="photonTree" )
	arguments.add_argument( "--order", action="store_true" )
	opts = arguments.parse_args()

	compareTrees( "ht" )
	#compareTrees( "(photons[0].matchedJetIndex>=0)*jets[photons[0].matchedJetIndex].pt+(photons[0].matchedJetIndex<0)*photons[0].pt", "GJets_400_inf_V02.28__31_tree.root" )
	#compareTrees( "photons[0]._ptJet" )
	#compareTrees( "photons[0].ptJet()")
	#compareTrees( "photons[0].pt+photons[0].chargedIso+photons[0].photonIso+photons[0].neutralIso")
	#compareTrees( "met" )
	#compareTrees( "ht" )
