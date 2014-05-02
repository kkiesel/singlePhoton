#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

def getFromEvent( event ):
	index = event.photons[0].matchedJetIndex
	return event.jets[index].pt
	#return event.jets[index].neutralEmEnergy
	#return event.jets[index].chargedHadronEnergy
	if index >=0:
		return event.jets[index].pt
	else:
		return event.photons[0].pt

def compareTrees( plot="photons.pt", filename="slimQCD_V02.28_tree.root" ):
	label, unit, binning = readAxisConf( "photons[0].ptJet()" )
	#binning = range(0,70,10) + binning
	#binning = range(0,1000, 30 )
	import array
	gH = ROOT.TH1F( randomName(), ";jet_{x};", len(binning)-1, array.array('d', binning) )
	fH = gH.Clone( randomName() )
	fH.SetLineColor(2)
	fH.SetMarkerColor(2)
	#cut = "!photons[0].isGen(0)"
	#cut = "1"

	gTree = readTree( filename, "photonTree" )
	fTree = readTree( filename, "photonJetTree" )

	#gH = getHisto( gTree, plot, color=1, cut=cut,firstBin=0,lastBin=1400 )
	#fH = getHisto( fTree, plot, color=2, cut=cut,firstBin=0,lastBin=1400 )

	for h, tree in [(gH, gTree), (fH, fTree)]:
		h.Sumw2()
		for event in tree:
			h.Fill( getFromEvent( event ), event.weight )


	for h in [gH, fH]:
		h.Scale( 1, "width" )

	mh = Multihisto()
	mh.addHisto( gH, "#gamma", draw="hist e" )
	mh.addHisto( fH, "#gamma_{jet}", draw="hist e" )

	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()


	from myRatio import Ratio
	r = Ratio( "#gamma/#gamma_{jet}", gH, fH )
	r.draw(None,2)


	SaveAs( can, "compare_test" )


def drawPhi( filenames ):
	t = readTree( filenames[0], "photonTree" )
	for f in filenames[1:]:
		t.Add( "%s/photonTree"%f )

	absolute = True

	label, unit, binning = readAxisConf( "photons[0].phi" )

	if absolute:
		binning = [ x for x in binning if x >= 0 ]

	import array
	gH1 = ROOT.TH1F( randomName(), ";#Delta#phi(#gamma,#slash{E}_{T});", len(binning)-1, array.array('d', binning) )
	gH2 = ROOT.TH1F( randomName(), ";#Delta#phi(#gamma,#slash{E}_{T});", len(binning)-1, array.array('d', binning) )
	gH3 = ROOT.TH1F( randomName(), ";#Delta#phi(#gamma,#slash{E}_{T});", len(binning)-1, array.array('d', binning) )
	gH2.SetLineColor(2)
	gH3.SetLineColor(3)
	for h in gH1, gH2, gH3:
		h.Sumw2()

	from math import fabs
	for e in t:

		dPhi = e.photons[0].DeltaPhi( e.metPhi )
		if absolute:
			dPhi = fabs( dPhi )

		if e.met<10:
			gH1.Fill( dPhi )
		elif e.met < 100:
			gH2.Fill( dPhi )
		else:
			gH3.Fill( dPhi )

	ft = readTree( filenames[0], "photonJetTree" )
	for f in filenames[1:]:
		ft.Add( "%s/photonJetTree"%f )

	cutStr = "photons[0].chargedIso/100<2.6 && photons[0].neutralIso/100<3.5+0.04*photons[0].pt && photons[0].photonIso/100<1.3+0.005*photons[0].pt"
	cutStr += "&& (photons[0].chargedIso>0 || photons[0].neutralIso>0 || photons[0].photonIso>0)"
	ft = ft.CopyTree( cutStr )

	fH1 = ROOT.TH1F( randomName(), ";#Delta#phi(#gamma,#slash{E}_{T});", len(binning)-1, array.array('d', binning) )
	fH2 = ROOT.TH1F( randomName(), ";#Delta#phi(#gamma,#slash{E}_{T});", len(binning)-1, array.array('d', binning) )
	fH3 = ROOT.TH1F( randomName(), ";#Delta#phi(#gamma,#slash{E}_{T});", len(binning)-1, array.array('d', binning) )
	fH2.SetLineColor(2)
	fH3.SetLineColor(3)
	for h in fH1, fH2, fH3:
		h.Sumw2()
		h.SetLineStyle(3)

	for e in ft:

		dPhi = e.photons[0].DeltaPhi( e.metPhi )
		if absolute:
			dPhi = fabs( dPhi )

		if e.met<10:
			fH1.Fill( dPhi )
		elif e.met < 100:
			fH2.Fill( dPhi )
		else:
			fH3.Fill( dPhi )


	for h in gH1, gH2, gH3, fH1, fH2, fH3:
		if h.Integral():
			h.Scale( 1./h.Integral(), "width" )

	mh = Multihisto()
	mh.addHisto( gH1, "met<10", draw="hist " )
	mh.addHisto( fH1, "met<10, loose", draw="hist " )
	mh.addHisto( gH2, "10<met<100", draw="hist " )
	mh.addHisto( fH2, "10<met<100, loose", draw="hist " )
	mh.addHisto( gH3, "100<met", draw="hist " )
	mh.addHisto( fH3, "100<met, loose", draw="hist " )


	c = ROOT.TCanvas()
	c.cd()
	mh.Draw()
	c.SaveAs("plots/deltaPhi_%s.pdf"%getSaveNameFromDatasets( filenames ))


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", default="met" )
	arguments.add_argument( "--tree", default="photonTree" )
	arguments.add_argument( "--order", action="store_true" )
	opts = arguments.parse_args()

	drawPhi( opts.filenames )

	#compareTrees( "photons[0].pt", "slimQCD_V02.28_tree.root" )
	#compareTrees( "photons[0].pt", "slimAllQCD_V02.28_tree.root" )
	#compareTrees( "photons[0].pt", "slimQCD_1000_inf_V02.28_tree.root" )
	#compareTrees( "(photons[0].matchedJetIndex>=0)*jets[photons[0].matchedJetIndex].pt+(photons[0].matchedJetIndex<0)*photons[0].pt", "GJets_400_inf_V02.28__31_tree.root" )
	#compareTrees( "photons[0]._ptJet" )
	#compareTrees( "photons[0].ptJet()")
	#compareTrees( "photons[0].pt+photons[0].chargedIso+photons[0].photonIso+photons[0].neutralIso")
	#compareTrees( "met" )
	#compareTrees( "ht" )
