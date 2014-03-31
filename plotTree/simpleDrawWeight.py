#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

st = ROOT.gROOT.GetStyle("tdrStyle")
paperWidth = 14.65 #cm


st.SetPaperSize(paperWidth/2,50.)

def getMixedWeigthHisto( filenames, predFilenames, commonCut, control=True, fillEmptyBins=False ):
	"""Calculate #photons/#photonFakes in bins of photons.ptJet and a second
	(global) variable.

	filenames: files containing photons
	predFilenames: files containing fakes
	"""

	regionCut = "1"

	xVar = "ht"
	yVar = "ht"
	xlabel, xunit, xbinning = "", "", range(500, 3001, 5 )
	#xlabel, xunit, xbinning = "", "", [0,500, 3000]
	ylabel, yunit, ybinning = "", "", [0, 3000]

	numerator = None
	for fileName in filenames:
		gTree = readTree( fileName, "photonTree" )
		num = createHistoFromTree2D( gTree, yVar+":"+xVar, "weight*( %s && %s )"%(regionCut, commonCut), xbinning, ybinning )
		if numerator:
			numerator.Add( num )
		else:
			numerator = num

	denominator = None
	for fileName in predFilenames:
		foTree = readTree( fileName, "photonJetTree" )
		den = createHistoFromTree2D( foTree, yVar+":"+xVar, "weight*( %s && %s )"%(regionCut, commonCut), xbinning, ybinning )
		if denominator:
			denominator.Add( den )
		else:
			denominator = den

	weight2D = divideHistos( numerator, denominator )

	# Set the weight and error for empty bins to one.
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if fillEmptyBins and not weight2D.GetBinContent( i, j ):
				weight2D.SetBinContent( i, j, 1 )
				weight2D.SetBinError( i, j, 1 )

	return weight2D


def writeWeight2DToFile( fileName, tree, h_weight, weightTreeName ):
	"""Write weight for a tree in another tree in a given file.
	This tree can be added to the original tree via 'AddFriend()'.

	fileName: name of file to which tree is written
	tree: tree which is weighted
	h_weight: two dimensional histogram with weights
	weighTreeName: name of the new tree
	"""
	weightTree = ROOT.TTree( weightTreeName, "Tree containing QCD weights" )
	import numpy
	weight = numpy.zeros( 1, dtype=float)
	weight_error = numpy.zeros( 1, dtype=float)
	weightTree.Branch( "w_compare", weight, "w_compare/D" )
	weightTree.Branch( "w_compare_error", weight_error, "w_compare_error/D" )

	from sys import stdout
	for event in tree:
		if not event.GetReadEntry()%10000:
			stdout.write( "\r%s / %s"%(event.GetReadEntry(), event.GetEntries() ) )
			stdout.flush()

		b = h_weight.FindBin( event.ht, event.ht )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()
	print

	f = ROOT.TFile( fileName, "update" )
	f.cd()
	weightTree.Write()
	f.Close()

def qcdPredictionHistos( filenames, plot, cut, modifyEmptyBins ):
	fHist, sysHist = None, None
	#cut+=" && !photons.isGen(2) && @photons.size()"
	for filename in filenames:
		fTree = readTree( filename, "photonJetTree" )
		fTree.AddFriend( "foWeightsCompare", filename )

		hist = getHisto( fTree, plot, weight="weight*w_compare", cut=cut, color=2, nBins=range(100,800,20)+range(800,1000,100), fillEmptyBins=modifyEmptyBins )

		sHist = getHisto( fTree, plot, weight="weight*w_compare_error", cut=cut, color=2, nBins=range(100,800,20)+range(800,1000,100), fillEmptyBins=modifyEmptyBins )


		if fHist:
			fHist.Add( hist )
			sysHist.Add( sHist )
		else:
			fHist = hist
			sysHist = sHist

	sysHist.SetFillColor( sysHist.GetLineColor() )
	sysHist.SetLineColor( sysHist.GetLineColor() )
	sysHist.SetFillStyle(3254)
	sysHist.SetMarkerSize(0)

	#fHist.SetLineColor(7)

	return fHist, sysHist



def getHistoFromFiles( plot, treename, filenames, cut, color=1 ):
	for filename in filenames:
		tree = readTree( filename, treename )
		h = getHisto( tree, plot, cut=cut, color=color, nBins=range(100,800,20)+range(800,1000,100) )
		#h = getHisto( tree, plot, cut=cut, color=color, nBins=range(500,1500,50) )
		if filename == filenames[0]:
			histo = h
		else:
			histo.Add( h )
	return histo


def compareTrees( plot="photons.pt", filenames=["slimAllQCD_V02.28_tree.root"], drawRatio=False ):
	cut = "1"
	#cut = plot+">110"
	#cut = "!photons[0].isGen(0)"
	#cut = "ht < 600"

	gH = getHistoFromFiles( plot, "photonTree", filenames, cut )
	fH = getHistoFromFiles( plot, "photonJetTree", filenames, cut, color=2 )

	weight2D = getMixedWeigthHisto( filenames, filenames, cut, control=True, fillEmptyBins=False )
	for filename in filenames:
		fTree = readTree( filename, "photonJetTree" )
		writeWeight2DToFile( filename, fTree, weight2D, "foWeightsCompare" )
	fgammaHist, fgammaWeightError = qcdPredictionHistos( filenames, plot, cut, True )
	for bin in range( fgammaHist.GetNbinsX()+2 ):
		fgammaHist.SetBinError( bin, sqrt( fgammaHist.GetBinError(bin)**2 + fgammaWeightError.GetBinContent(bin)**2 ) )

	fH = fgammaHist


	for h in [gH, fH]:
		h.Scale( 1./h.Integral() )
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


	mh = Multihisto()
	mh.leg.SetX1(0.608)
	mh.leg.SetX2(0.95)
	mh.addHisto( gH, "#gamma_{#text{tight}}", draw="hist e" )
	mh.addHisto( fH, "#gamma_{#text{loose}}", draw="hist e" )

	can = ROOT.TCanvas()
	#can.SetBottomMargin(0)
	#can.SetTopMargin(0)
	#can.SetRightMargin(0)
	can.SetLeftMargin(0.19)
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
	SavePad( "compare2_%s_%s"%(plot,shortName( filenames )) )
	ROOT.gPad.SaveAs("~/master/documents/thesis/plots/compare_%s_%s.tex"%(manipulateSaveName(plot),shortName( filenames )) )
	correctTiksPlot("/home/knut/master/documents/thesis/plots/compare_%s_%s.tex"%(manipulateSaveName(plot),shortName( filenames )) )
	ROOT.SetOwnership( can, False )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	#arguments.add_argument( "--plots", nargs="+", default=["photons[0].pt"] )
	arguments.add_argument( "--plots", nargs="+", default=["photons[0].ptJet()", "photons[0].pt", "ht"] )
	arguments.add_argument( "--ratio", action="store_true" )
	opts = arguments.parse_args()

	for plot in opts.plots:
		compareTrees( plot, opts.filenames, opts.ratio )

