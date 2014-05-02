#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from sys import stdout

# variable to reweight
reweightVar = "ht"

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
	weightTree.Branch( "w_qcd", weight, "w_qcd/D" )
	weightTree.Branch( "w_qcd_error", weight_error, "w_qcd_error/D" )

	for event in tree:
		if not event.GetReadEntry()%10000:
			stdout.write( "\r%s / %s"%(event.GetReadEntry(), event.GetEntries() ) )
			stdout.flush()

		b = h_weight.FindBin( event.photons.at(0).ptJet(), eval("event.%s"%reweightVar) )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()
	print

	f = ROOT.TFile( fileName, "update" )
	f.cd()
	weightTree.Write()
	f.Close()

def getWeightHisto2D( gControlTree, foControlTree, commonCut, datasetAbbr ):
	"""The histogram for the weights in created here."""

	xVar = "photons[0].ptJet()"
	yVar = reweightVar
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

	weight_numerator = createHistoFromTree2D( gControlTree, yVar+":"+xVar, "weight*(met<100&&%s)"%commonCut, xbinning, ybinning )
	weight_denominator = createHistoFromTree2D( foControlTree, yVar+":"+xVar, "weight*(met<100&&%s)"%commonCut, xbinning, ybinning )
	weight2D = divideHistos( weight_numerator, weight_denominator )

	# Set the weight and error for empty bins to one.
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if not weight2D.GetBinContent( i, j ):
				weight2D.SetBinContent( i, j, 1 )
				weight2D.SetBinError( i, j, 1 )

	return weight2D


def ePred( filename, color, cut="photons.isGen(0) && @photons.size()" ):
	tree = readTree( filename, "photonElectronTree" )
	h = getHisto( tree, "met", cut=cut, color=color )
	return h

def finalDistribution( data, tt, wj, gw, gz ):
	gDataTree = readTree( data, "photonTree" )
	eDataTree = readTree( data, "photonElectronTree" )
	fDataTree = readTree( data, "photonJetTree" )

	signal = getHisto( gDataTree, "met", color=1 )
	eHisto = applyFakeRateEWK( getHisto( eDataTree, "met", color=8 ), 0.015 )

	weights = getWeightHisto2D( gDataTree, fDataTree, "1", "" )
	writeWeight2DToFile( data, fDataTree, weights, "foWeights" )
	fDataTree.AddFriend( "foWeights", data )
	fHisto = getHisto( fDataTree, "met", weight="weight*w_qcd", color=38 )
	fHistoUp = getHisto( fDataTree, "met", weight="weight*(w_qcd+w_qcd_error)", color=38 )
	fHistoDown = getHisto( fDataTree, "met", weight="weight*(w_qcd-w_qcd_error)", color=38 )
	sysHist = fHistoUp.Clone( randomName() )
	for bin in range( fHistoUp.GetNbinsX()+2 ):
		up = fHistoUp.GetBinContent(bin)
		down = fHistoDown.GetBinContent(bin)
		sysHist.SetBinContent( bin, (up+down)/2 )
		sysHist.SetBinError( bin, (up-down)/2 )
		fHisto.SetBinError( bin, sqrt( fHisto.GetBinError(bin)**2 + sysHist.GetBinError(bin)**2 ) )


	mh = Multihisto()
	mh.addHisto( fHisto, "QCD(+#gamma)", True )
	mh.addHisto( eHisto, "e#rightarrow#gamma", True )
	#mh.addHisto( ePred( tt, ROOT.kOrange ), "#gammat#bar{t}", True )
	#mh.addHisto( ePred( wj, ROOT.kOrange-2 ), "#gammaWJet", True )
	#mh.addHisto( ePred( gz, ROOT.kOrange+1 ), "#gammaZ", True )
	#mh.addHisto( ePred( gw, ROOT.kOrange+2 ), "#gammaW", True )
	mh.addHisto( signal, "Data", draw="pe0" )

	can = ROOT.TCanvas()
	hPad = ROOT.TPad("hPad", "Histogram", 0, 0.2, 1, 1)
	hPad.cd()
	mh.Draw()
	info = PlotCaption()
	info.Draw()

	ratioPad = ROOT.TPad("ratioPad", "Ratio", 0, 0, 1, 0.2)
	ratioPad.cd()
	ratioPad.SetLogy( False )
	from myRatio import Ratio
	r = Ratio( "Data/Bkg", signal, mh.histos[-1][0].GetStack().Last() )
	ratio, sys, one = r.draw(.5,1.5)
	ratio.Draw("same e0")
	sys.Draw("same e2")
	one.Draw()

	can.cd()
	hPad.Draw()
	ratioPad.Draw()
	SaveAs( can, "test" )
	ROOT.SetOwnership( hPad, False )
	ROOT.SetOwnership( ratioPad, False )



if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	#arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--plot", nargs="+", default = ["met"] )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "photons[0].ptJet()","Length$(jets.pt)", "Length$(photons.pt)"]

	version = "19"
	data = "PhotonHad_V02.%s_tree.root"%version
	tt = "slimTTJets_V02.%s_tree.root"%version
	wj = "slimWJets_V02.%s_tree.root"%version
	gw = "slimWGamma_V02.%s_tree.root"%version
	gz = "slimZGamma_V02.%s_tree.root"%version

	finalDistribution( data, tt, wj, gw, gz )

