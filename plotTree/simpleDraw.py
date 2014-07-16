#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

st = ROOT.gROOT.GetStyle("tdrStyle")
paperWidth = 14.65 #cm
#st.SetPaperSize(paperWidth/2,50.)

st.SetPaperSize(paperWidth/2,50.)

def getHistoFromFiles( plot, treename, filenames, cut, color=1 ):
	for filename in filenames:
		tree = readTree( filename, treename )
		#h = getHisto( tree, plot, cut=cut, color=color, nBins=range(100,800,20)+range(800,1000,100) )
		h = getHisto( tree, plot, cut=cut, color=color )
		if filename == filenames[0]:
			histo = h
		else:
			histo.Add( h )
	return histo


def compareTrees( plot="met", filenames=["slimAllQCD_V02.28_tree.root"], drawRatio=False ):
	#cut = "@electrons.size()==0 && @muons.size()==0"
	#cut = plot+">110"
	#cut = "!photons[0].isGen(0)"
	#cut = "ht < 600"
	cut = "1"

	gH = getHistoFromFiles( plot, "photonTree", filenames, "1" )
	fH = getHistoFromFiles( plot, "photonJetTree", filenames, "1", color=2 )

	for h in [gH, fH]:
		#h.Scale( 1./h.Integral() )
		h.SetMarkerSize(0)
		h.GetXaxis().SetTitleOffset( 1.03 )
		h.GetYaxis().SetTitleOffset( 1.5 )
		h.GetXaxis().SetLabelSize( 0.0633790992496425 )
		h.GetXaxis().SetTitleSize( 0.0633790992496425 )
		h.GetYaxis().SetLabelSize( 0.0633790992496425 )
		h.GetYaxis().SetTitleSize( 0.0633790992496425 )
		h.GetYaxis().SetLabelOffset(0)

		if plot == "photons[0].pt":
			h.GetXaxis().SetTitle( "$p_{T}$ [GeV]" )
		if plot == "photons[0].ptJet()":
			h.GetXaxis().SetTitle( "$p_{T^*}$ [GeV]" )

	ratio = gH.Clone( randomName() )
	ratio.Divide( fH )


	mh = Multihisto()
	#mh.leg.SetX1(0.608)
	#mh.leg.SetX2(0.95)
	mh.addHisto( gH, "#gamma_{#text{tight}}", draw="hist e" )
	mh.addHisto( fH, "#gamma_{#text{loose}}", draw="hist e" )

	can = ROOT.TCanvas()
	#can.SetBottomMargin(0)
	#can.SetTopMargin(0)
	#can.SetRightMargin(0)
	#can.SetLeftMargin(0.19)
	can.cd()
	mh.Draw()

	pc1 = ROOT.TLatex(0,.96, "CMS Private Work")
	pc2 = ROOT.TLatex( .51,.96, "\SI{19.8}{fb^{-1}} #sqrt{s}=\SI{8}{TeV}")
	for pc in [pc1, pc2]:
		pc.SetNDC()
		pc.SetTextSize(0.06311227345609463)
		pc.Draw()


	if drawRatio:
		from myRatio import Ratio
		r = Ratio( "#gamma/#gamma_{jet}", gH, fH )
		r.draw()

	#can.SetFillColor(ROOT.kGreen)
	SavePad( "compare2_%s_%s"%("test",shortName( filenames )) )
	#ROOT.gPad.SaveAs("~/master/documents/thesis/plots/compare_%s_%s.tex"%(manipulateSaveName(plot),shortName( filenames )) )
	ROOT.SetOwnership( can, False )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	#arguments.add_argument( "--plots", nargs="+", default=["photons[0].pt"] )
	arguments.add_argument( "--plots", nargs="+", default=[ "met"] )
	arguments.add_argument( "--ratio", action="store_true" )
	opts = arguments.parse_args()

	for plot in opts.plots:
		compareTrees( plot, opts.filenames, opts.ratio )

