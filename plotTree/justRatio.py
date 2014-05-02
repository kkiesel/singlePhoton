#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

st = ROOT.gROOT.GetStyle("tdrStyle")


def texStyle():
	#style = ROOT.TStyle("texStyle", "Style to generate PGF/TikZ")
	style = ROOT.gStyle

	# Canvas
	style.SetCanvasColor( ROOT.kWhite )
	style.SetCanvasBorderMode(0)
	style.SetCanvasDefH(1000)
	style.SetCanvasDefW(2000)

	# Pad:
	style.SetPadBorderMode(0)

	# Margins:
	style.SetPadTopMargin(0.06)
	style.SetPadBottomMargin(0.13)
	style.SetPadLeftMargin(0.16)
	style.SetPadRightMargin(0.02)

	# axis titles
	style.SetTitleFont(42, "xyz")
	style.SetTitleSize(0.09, "xyz")
	style.SetLabelSize(0.09, "xyz")

	# Statistic box
	style.SetOptStat(0)

	# Other
	style.SetPalette(1)
	# x-value is approximation of text width, y value large
	style.SetPaperSize(14.65,50.)

	style.SetPadTickX(1)
	style.SetPadTickY(1)

	style.cd()
	return style

def getHistoFromFiles( plot, treename, filenames, cut, color=1 ):
	for filename in filenames:
		tree = readTree( filename, treename )
		try:
			tree.AddFriend( treename+"angular", filename )
		except:
			print "tree %s not avaiable in %s"%(treename+"angular", filename)

		if plot == "met":
			h = getHisto( tree, plot, cut=cut, color=color,nBins = range(0,101,10) )
		elif "Sin" in plot:
			h = getHisto( tree, plot, cut=cut, color=color,nBins = range(0,101,10) )
		elif "Cos" in plot:
			h = getHisto( tree, plot, cut=cut, color=color,nBins = range(-100,101,10) )
		else:
			h = getHisto( tree, plot, cut=cut, color=color )
		if filename == filenames[0]:
			histo = h
		else:
			histo.Add( h )
	return histo


def compareTrees( plot="met", filenames=["slimAllQCD_V02.28_tree.root"], drawRatio=False ):
	cut = "1"
	fcut = cut
	fcut += "&& photons[0].chargedIso/10<2.6 && photons[0].neutralIso/10<3.5+0.04*photons[0].pt && photonIso/10<1.3+0.005*photons[0].pt"
	fcut += "&& photons[0].chargedIso*10>2.6 && photons[0].neutralIso*10>3.5+0.04*photons[0].pt && photonIso*10>1.3+0.005*photons[0].pt"


	gH = getHistoFromFiles( plot, "photonTree", filenames, cut )
	fH = getHistoFromFiles( plot, "photonJetTree", filenames, fcut, color=1 )
	#fH = getHistoFromFiles( plot, "photonTree", filenames, cut+"&& fabs(photons[0].eta)>=%s"%etaCut )



	ratio = gH.Clone( randomName() )
	ratio.Divide( fH )
	ratio.SetYTitle( "#gamma_{tight}/#gamma_{loose}" )


	can = ROOT.TCanvas()
	can.SetLogy(0)
	can.cd()
	ratio.Draw(" e")

	pc1 = ROOT.TLatex(0,.96, "CMS Private Work")
	pc2 = ROOT.TLatex( .42,.96, "\SI{19.8}{fb^{-1}} #sqrt{s}=\SI{8}{TeV}#, #geq1#ggamma(#geq1#fgamma),#geq2#text{jets}")
	for pc in [pc1, pc2]:
		pc.SetNDC()
		pc.SetTextSize(0.06311227345609463)
		pc.Draw()


	if drawRatio:
		from myRatio import Ratio
		r = Ratio( "#gamma/#gamma_{jet}", gH, fH )
		r.draw()

	#can.SetFillColor(ROOT.kGreen)
	SavePad( "ratioCompare2_%s_%s"%(plot,shortName( filenames )) )
	#ROOT.gPad.SaveAs("~/master/documents/thesis/plots/ratioCompare_%s_%s.tex"%(manipulateSaveName(plot),shortName( filenames )) )
	#ROOT.SetOwnership( can, False )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	#arguments.add_argument( "--plots", nargs="+", default=["photons[0].pt"] )
	arguments.add_argument( "--plots", nargs="+", default=[ "met"] )
	arguments.add_argument( "--ratio", action="store_true" )
	opts = arguments.parse_args()

	for plot in opts.plots:
		compareTrees( plot, opts.filenames, opts.ratio )

