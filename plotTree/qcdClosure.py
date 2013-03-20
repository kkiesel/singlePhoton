#! /usr/bin/env python2
# -*- coding: utf-8 -*-

import ROOT
import argparse
import ConfigParser
import os.path
import Styles
Styles.tdrStyle()

axisConf = ConfigParser.SafeConfigParser()
axisConf.read("axis.cfg")

integratedLuminosity = 19.3 #fb

def readTree( filename, treename = "photonTree" ):
	"""
	filename: name of file containing the tree
	treename: name of the tree
	returns: TTree Object
	"""
	if not os.path.isfile(filename):
		print( "File %s does not exist"%filename)
	tree = ROOT.TChain( treename )
	tree.AddFile( filename )
	return tree

def createHistoFromTree2D(tree, variable, weight="", nBin=50):
	"""
	tree: tree to create histo from
	variable: variable to plot (must be a branch of the tree)
	weight: weights to apply (e.g. "var1*(var2 > 15)" will use weights from var1 and cut on var2 > 15
	nBins, firstBin, lastBin: number of bins, first bin and last bin (same as in TH1F constructor)
	nEvents: number of events to process (-1 = all)
	returns: histogram
	"""
	from random import randint
	from sys import maxint
	#make a random name you could give something meaningfull here,
	#but that would make this less readable
	name = "%x"%(randint(0, maxint))

	# automatic binning if lastbin < firstbin
	result = ROOT.TH2F(name, variable, nBin, 0, -1, nBin, 0, -1)
	result.Sumw2()
	tree.Draw("%s>>%s"%(variable, name), weight, "goff,colz")
	return result

def createHistoFromTree(tree, variable, weight="", nBins=100, firstBin=None, lastBin=None, nEvents=-1):
	"""
	tree: tree to create histo from
	variable: variable to plot (must be a branch of the tree)
	weight: weights to apply (e.g. "var1*(var2 > 15)" will use weights from var1 and cut on var2 > 15
	nBins, firstBin, lastBin: number of bins, first bin and last bin (same as in TH1F constructor)
	nEvents: number of events to process (-1 = all)
	returns: histogram
	"""
	if ":" in variable:
		return createHistoFromTree2D(tree, variable, weight="")
	from ROOT import TH1F
	from random import randint
	from sys import maxint
	if nEvents < 0:
		nEvents = maxint
	if firstBin == None:
		firstBin = tree.GetMinimum( variable )
	if lastBin == None:
		lastBin = tree.GetMaximum( variable )
	#make a random name you could give something meaningfull here,
	#but that would make this less readable
	name = "%x"%(randint(0, maxint))
	result = TH1F(name, variable, nBins, firstBin, lastBin)
	result.Sumw2()
	tree.Draw("%s>>%s"%(variable, name), weight, "goff", nEvents)
	return result

def plot2D( tree, plot, cut, save=False ):
	# read axis config file
	var = plot.split(":")
	try:
		xlabel = axisConf.get( var[0], "label" )
	except:
		xlabel = ""
	try:
		ylabel = axisConf.get( var[1], "label" )
	except:
		ylabel = ""
	try:
		xunit = axisConf.get( var[0], "unit" )
	except:
		xunit = ""
	try:
		yunit = axisConf.get( var[1], "unit" )
	except:
		yunit = ""
	if xunit != "":
		xunit = " ["+xunit+"]"
	if yunit != "":
		yunit = " ["+yunit+"]"

	# modify histo
	histo = createHistoFromTree2D( tree, plot, cut )
	histo.SetTitle(";%s%s;%s%s"%(ylabel,yunit,xlabel,xunit))

	# draw histo
	canvas = ROOT.TCanvas()
	canvas.cd()
	histo.Draw("colz")

	# draw information
	cutText = ROOT.TPaveText(.5,.8,.9,.9,"ndc")
	cutText.SetBorderSize(0)
	cutText.SetFillColor(0)
	for line in cut.split("&&"):
		cutText.AddText( line )
	cutText.Draw()

	topInfo = ROOT.TPaveText(.1,.95,1,1,"ndc")
	topInfo.SetBorderSize(0)
	topInfo.SetFillColor(0)
	topInfo.AddText("Work in progress, %sfb^{-1}, #sqrt{s}=8TeV"%"?")
	topInfo.Draw()

	if save:
		canvas.SaveAs("%.pdf"%variable)
	else:
		raw_input()

def fillWeights( tree, jetCut, photonCut ):

	# this is the binning with witch the fo reweighting will be done and has to be studied
	nBins = 50
	xMin = 80
	xMax = 300
	h_photonPt = createHistoFromTree( tree, "photon.pt", photonCut, nBins, xMin, xMax )
	h_jetPt = createHistoFromTree( tree, "photon.pt", jetCut, nBins, xMin, xMax )

	# h_jetPt is now histogram with w^{-1} = h_jetPt / h_photonPt
	h_jetPt.Divide( h_photonPt )

	for event in range( tree.GetNEntries() ):
		tree.GetEntry( event )
		tree.weight *= h_jetPt.GetBinContent( h_jetPt.FindBin( event.photon.pt ) )



def plot( tree, plot, cut, save=False ):
	# read axis config file
	try:
		label = axisConf.get( plot, "label" )
		unit = axisConf.get( plot, "unit" )
		if unit != "":
			unit = " ["+unit+"]"
	except:
		print "please specify label and unit in axis configuration file"
		label = plot
		unit = ""

	# modify histo
	histo = createHistoFromTree( tree, plot, cut )
	histo.SetTitle(";%s%s;Entries"%(label,unit))
	histo.SetLineColor(1)

	# draw histo
	canvas = ROOT.TCanvas()
	canvas.cd()
	canvas.SetLogy()
	histo.Draw("hist")

	# draw information
	cutText = ROOT.TPaveText(.5,.8,.9,.9,"ndc")
	cutText.SetBorderSize(0)
	cutText.SetFillColor(0)
	for line in cut.split("&&"):
		cutText.AddText( line )
	cutText.Draw()

	topInfo = ROOT.TPaveText(.1,.95,1,1,"ndc")
	topInfo.SetBorderSize(0)
	topInfo.SetFillColor(0)
	topInfo.AddText("Work in progress, %sfb^{-1}, #sqrt{s}=8TeV"%"?")
	topInfo.Draw()

	if save:
		canvas.SaveAs("%.pdf"%variable)
	else:
		raw_input()


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("-f", "--file", default="myTree.root", help="ROOT file containing tree")
	arguments.add_argument("-c", "--cut", default="", help="Cut string")
	arguments.add_argument("-d", "--distribution", default =['met'], nargs="+",
			help="Distribution which shall be plotted.")
	arguments.add_argument("--save", action="store_true", help="Save canvas as pdf.")

	opts = arguments.parse_args()

	tree = readTree( opts.file )

	# see https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonID2012
	# for more information on 2012 Photon ID
	# TODO: Conversion safe electron veto
	# TODO: Single tower H/E
	# TODO: correct isolation with ρ
	photonCut2012ID = "photon.sigmaIetaIeta<0.012 \
&& photon.chargedIso < 2.6 \
&& photon.neutralIso < 3.5 + 0.04*photon.pt \
&& photon.photonIso < 1.3 + 0.005*photon.pt"

	jetPhotonCut = "photon.sigmaIetaIeta > 0.012 && photon.sigmaIetaIeta<0.02 \
&& photon.chargedIso < 2.6 \
&& photon.neutralIso < 3.5 + 0.04*photon.pt \
&& photon.photonIso < 1.3 + 0.005*photon.pt"

	fillWeights( tree, jetPhotonCut, photonCut2012ID )

	if opts.cut == "photon":
		opts.cut = photonCut2012ID
	elif opts.cut == "jetphoton":
		opts.cut = jetPhotonCut

	if "all" in opts.distribution:
		opts.distribution = axisConf.sections()

	for plots in opts.distribution:
		if ":" in plots:
			plot2D( tree, plots, opts.cut, opts.save )
		else:
			plot( tree, plots, opts.cut, opts.save )

