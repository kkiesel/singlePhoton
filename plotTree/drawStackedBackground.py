#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
ROOT.gStyle.SetPadRightMargin(0.05)

colors = {
		"GJets": ROOT.kCyan-7,
		"GJets_200_400": ROOT.kCyan-10,
		"GJets_400_inf":ROOT.kCyan-7,

		"QCD": ROOT.kCyan+3,
		"QCD_250_500": ROOT.kCyan+3,
		"QCD_500_1000": ROOT.kCyan+2,
		"QCD_1000_inf": ROOT.kCyan+1,

		"WJets": ROOT.kGreen,
		"WJets_250_inf": ROOT.kGreen,
		"WJets_250_300": ROOT.kGreen,
		"WJets_300_400": ROOT.kGreen+1,
		"WJets_400_inf": ROOT.kGreen+2,

		"TTJets": ROOT.kGreen-3,
		"TTFull": ROOT.kGreen-2,
		"TTSemi": ROOT.kGreen-1,
		"TTHadronic": ROOT.kGreen-0,
		"WGamma": ROOT.kRed-7,
		"WGamma_50_130": ROOT.kRed-7,
		"WGamma_130_inf": ROOT.kRed-6,
		"WGamma_50_inf": ROOT.kRed-6,
		"ZGamma": ROOT.kRed-4,
		"ZGammaNuNu": ROOT.kRed-3,
		"ZGammaLL": ROOT.kRed-2,
		"TTGamma": ROOT.kRed,
		}

combinedDatasets = {
		"GJets": ["GJets_200_400", "GJets_400_inf"],
		"QCD": ["QCD_250_500", "QCD_500_1000", "QCD_1000_inf"],
		"W": ["WJets_250_300", "WJets_300_400", "WJets_400_inf"],
		"WGamma": ["WGamma_50_130", "WGamma_130_inf"]
	}

def drawStackedBackground( plot, treeName, listOfFiles, sumBinned, order=False ):

	cut = "!@electrons.size() && !@muons.size()"

	# if a data histogram is present, the mc integral will scaled to data
	mcHists = []
	dataHists = []
	for fileName in listOfFiles:
		datasetAbbr = getDatasetAbbr( fileName )
		try:
			color = colors[datasetAbbr]
		except:
			color = 1
		tree = readTree( fileName, treeName )
		histo = getHisto( tree, plot, color=color, cut=cut )
		if "PhotonHad" in fileName:
			dataHists.append( (datasetAbbr, histo ) )
		else:
			mcHists.append( ( datasetAbbr, histo ) )

	# merge datahists if there
	dataHist = addHistos( [histo for datasetAbbr, histo in dataHists] ) if dataHists else None

	scale = 1.*dataHist.Integral()/addHistos( [histo for datasetAbbr, histo in mcHists] ).Integral() if dataHists else 1
	scale = 1.

	# combine MC datasets
	if sumBinned:
		for combiAbbr, abbrList in combinedDatasets.iteritems():
			if set(abbrList).issubset( set( [ a for a, h in mcHists ] ) ):
				indices = []
				for a,b in mcHists:
					if a in abbrList:
						indices.append(mcHists.index((a,b)))
				#print "a valid combination was found"
				histosToAdd = [ h for a,h in mcHists if a in abbrList ]
				thisSum = addHistos( histosToAdd )
				mcHists = [ (a,h) for a,h in mcHists if a not in abbrList]
				mcHists.insert(min(indices),(combiAbbr,thisSum))

	mh = Multihisto()
	mh.setMinimum(0.01)
	mh.leg.SetX1(0.5)
	mh.leg.SetX2(0.95)
	mh.leg.SetY1(0.65)
	mh.leg.SetFillStyle(0)
	mh.leg.SetNColumns(2)
	mh.leg.SetColumnSeparation(-0.50)
	mh.orderByIntegral = order
	if dataHist:
		mh.addHisto( dataHist, "Data", draw="ep" )

	for abbr, hist in mcHists:
		hist.Scale( scale )

		hist.GetXaxis().SetTitle("#met#text{ [GeV]}")
		hist.GetXaxis().SetNdivisions(5,5,0, False)
		for a in hist.GetXaxis(), hist.GetYaxis():
			a.SetLabelSize(1./31.4485/0.425531/0.669117)
			a.SetTitleSize(1./31.4485/0.425531/0.80294)
			a.SetLabelOffset(1000)
		hist.GetXaxis().SetTitleOffset(20.0)
		hist.GetYaxis().SetTitleOffset(1.6)

		mh.addHisto( hist, datasetToLatex(abbr), toStack=True, draw="hist" )

	#totalBG = addHistos( [ h for abbr, h in mcHists ] )
	#totalBG.SetLineColor(2)
	#mh.addHisto( totalBG, "SM Simulation", toStack=False, draw="hist" )

	# add signal histos
	#signalFiles = ["slimW_1000_1020_375_V02.44_tree.root", "slimW_1200_1120_375_V02.44_tree.root"]
	signalFiles = ["slimW_1200_1120_375_V02.44_tree.root"]
	for iColor, sf in enumerate(signalFiles):
		tree = readTree( sf, treeName )
		histo = getHisto( tree, plot, color=ROOT.kBlue+iColor, cut=cut )
		#mh.addHisto( histo, datasetToLatex(getDatasetAbbr(sf)), draw="hist" )
		#mh.addHisto( histo, "Signal", draw="hist" )


	ROOT.gStyle.SetPaperSize(14.6/2.35,50.)
	ROOT.gStyle.SetPadTopMargin(0.05)
	ROOT.gStyle.SetPadRightMargin(0.02)
	ROOT.gStyle.SetPadLeftMargin(0.17)

	infoText = ROOT.TLatex(0,.96, "#text{CMS Private Work  }#geq1#gamma_{#text{pixel}},#geq2#text{jets}" )
	infoText.SetNDC()
	infoText.SetTextSize(1./14.3823)

	can = ROOT.TCanvas()
	can.cd()
	mh.Draw()
	xax = mh.stack.GetXaxis()
	xax.SetNdivisions(5,5,0,False)
	for a in xax, mh.stack.GetYaxis():
		a.SetTitleSize(1./13.3823)
		a.SetLabelSize(1./13.3823)

	infoText.Draw()
	if dataHist:
		from myRatio import Ratio
		den = totalBG
		#den = mh.stack.GetStack().Last().Clone( randomName() )
		den.SetLineColor(2)
		r = Ratio( "Data/Sim.", dataHist, den )
		r.draw(0,2)

	allDatasetAbbr = getSaveNameFromDatasets( listOfFiles )
	SaveAs( can, "stackedHisto_%s_%s_%s"%(treeName, plot,allDatasetAbbr) )
	can.SaveAs("~/master/documents/thesis/plots/stackedHistos_%s_%s_%s.tex"%(treeName,plot,allDatasetAbbr) )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--plot", nargs="+", default=["met"] )
	arguments.add_argument( "--tree", nargs="+", default=["photonTree"] )
	arguments.add_argument( "--order", action="store_true" )
	arguments.add_argument( "--sum", action="store_true" )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "nGoodJets", "photons[0].ptJet()", "photons[0].pt", "photons[0].eta", "photons[0].phi", "photons[0].sigmaIetaIeta", "photons[0].chargedIso", "photons[0].neutralIso", "photons[0].photonIso"]

	if opts.tree == ["all"]:
		opts.tree = [ "photonTree", "photonJetTree", "photonElectronTree" ]

	for plot in opts.plot:
		for tree in opts.tree:
			drawStackedBackground( plot, tree, opts.filenames, opts.sum, opts.order )
