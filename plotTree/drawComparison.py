#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

class Dataset:
	def __init__( self, listOfFileNames, color ):
		self.listOfFileNames = listOfFileNames
		self.datasetAbbr = mergeDatasetAbbr( [ getDatasetAbbr(x) for x in listOfFileNames ] )[0]
		self.color = color

	def plot( self, treeName="photonTree", plot="met", cut="photons[0].ptJet()>=100",color=None ):
		hist = None
		if color:
			self.color = color
		for filename in self.listOfFileNames:
			tree = readTree( filename, treeName )
			h = getHisto( tree, plot, cut=cut, color=self.color )
			if hist:
				hist.Add( h )
			else:
				hist = h
		return hist

gjets = Dataset( ["slimGJets_200_400_V02.35_tree.root", "slimGJets_400_inf_V02.35_tree.root"], ROOT.kCyan-7 )
qcd = Dataset( ["slimQCD_250_500_V02.35_tree.root", "slimQCD_500_1000_V02.35_tree.root", "slimQCD_1000_inf_V02.35_tree.root" ], ROOT.kCyan+3 )
tt = Dataset( ["slimTTJets_V02.35_tree.root"], ROOT.kGreen-3 )
wjets = Dataset( ["slimWJets_250_300_V02.35_tree.root","slimWJets_300_400_V02.35_tree.root", "slimWJets_400_inf_V02.35_tree.root"], ROOT.kGreen )
gw = Dataset( ["slimWGamma_50_130_V02.35_tree.root","slimWGamma_130_inf_V02.35_tree.root"], ROOT.kRed )
gz = Dataset( ["slimZGamma_V02.35_tree.root"], ROOT.kRed-4 )


def fillMh( mh, datasetsToStack, treeName ):
	for dset in datasetsToStack:
		mh.addHisto( dset.plot(treeName), datasetToLatex(dset.datasetAbbr), toStack=True )
	return mh

def drawMCStack( treeName="photonTree" ):
	mh = Multihisto()
	mh.orderByIntegral = False

	datasetsToStack = [gz,tt,wjets,qcd,gjets]

	mh = fillMh( mh, datasetsToStack, treeName )

	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()
	info = PlotCaption(treeName=treeName)
	info.Draw()
	allDatasetAbbr = ''.join([x.datasetAbbr for x in datasetsToStack ])
	plot = "met"
	SaveAs( can, "stackedHisto_%s_%s_%s"%(treeName, plot,allDatasetAbbr) )


drawMCStack()
drawMCStack("photonJetTree")
drawMCStack("photonElectronTree")

def drawGenComposition( treeName, sample ):
	mh = Multihisto()
	mh.addHisto( sample.plot(treeName, cut="photons[0].isGen(0)",color=2), "gen #gamma" )
	mh.addHisto( sample.plot(treeName, cut="photons[0].isGen(1)",color=3), "gen e" )
	mh.addHisto( sample.plot(treeName, cut="photons[0].isGen(2)",color=4), "gen hadron" )
	mh.addHisto( sample.plot(treeName, cut="1",color=1), "inclusive" )
	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()
	info = PlotCaption(treeName=treeName)
	info.Draw()
	allDatasetAbbr = sample.datasetAbbr
	plot = "met"
	SaveAs( can, "composition_%s_%s_%s"%(treeName, plot,allDatasetAbbr) )

def drawGenComposition2Samples( treeName, sample1, sample2 ):

	histos = {
			"W": sample1.plot(treeName, cut="1"),
			"W (gen #gamma)": sample1.plot(treeName, cut="photons[0].isGen(0)"),
			"W (gen e)": sample1.plot(treeName, cut="photons[0].isGen(1)"),
			"W (gen had)": sample1.plot(treeName, cut="photons[0].isGen(3)"),
			"#gammaW": sample2.plot(treeName, cut="1"),
			"#gammaW (gen #gamma)": sample2.plot(treeName, cut="photons[0].isGen(0)"),
			"#gammaW (gen e)": sample2.plot(treeName, cut="photons[0].isGen(1)"),
			"#gammaW (gen had)": sample2.plot(treeName, cut="photons[0].isGen(2)"),
			}

	mh = Multihisto()
	for name, h in histos.iteritems():
		if "(gen #gamma)" in name:
			h.SetLineStyle(2)
		if "(gen e)" in name:
			h.SetLineStyle(3)
		if "(gen had)" in name:
			h.SetLineStyle(4)

		mh.addHisto( h, name )
	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()
	info = PlotCaption(treeName=treeName)
	info.Draw()
	allDatasetAbbr = sample1.datasetAbbr+sample2.datasetAbbr
	plot = "met"
	SaveAs( can, "composition2samples_%s_%s_%s"%(treeName, plot,allDatasetAbbr) )


drawGenComposition( "photonTree", tt )
drawGenComposition( "photonTree", wjets )

drawGenComposition( "photonJetTree", tt )
drawGenComposition( "photonJetTree", wjets )


drawGenComposition2Samples( "photonTree", wjets, gw )
drawGenComposition2Samples( "photonJetTree", wjets, gw )


