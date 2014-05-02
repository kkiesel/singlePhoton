#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

st = ROOT.gROOT.GetStyle("tdrStyle")

def drange(start, stop, step):
	r = start
	while r < stop:
		yield r
		r += step

def stepRange( start, stop, nSteps=20):
	step = 1.*(stop-start)/nSteps
	return drange( start, stop, step )


def getHistoFromFiles( plot, treename, filenames, cut):
	for filename in filenames:
		tree = readTree( filename, treename )

		if filename.startswith("PhotonHad"):
			weight="1"
		else:
			weight="weight"

		if plot == "met" and True:
			h = getHisto( tree, plot, cut, weight=weight, nBins = range(0,101,10), appendOverflowBin=False )
		else:
			h = getHisto( tree, plot, cut, weight=weight, appendOverflowBin=False )
		if filename == filenames[0]:
			histo = h
		else:
			histo.Add( h )
	return histo

def compareHists( filenames, tight, fcut, draw=False ):
	loose = getHistoFromFiles( "met", "photonJetTree", filenames, fcut )
	loose.SetLineColor(2)
	if loose.Integral(0,loose.FindBin(99)):
		loose.Scale( tight.Integral(0,loose.FindBin(99)) / loose.Integral(0,loose.FindBin(99)) )
	option = "uu norm" if filenames[0].startswith("PhotonHad") else "ww"
	if draw:

		c = ROOT.TCanvas()
		c.cd()
		mh = Multihisto()
		mh.addHisto( tight.Clone( randomName() ), "#gamma_{tight}", draw="hist e")
		mh.addHisto( loose, "#gamma_{loose}", draw="hist e")
		mh.Draw()
		from myRatio import Ratio
		r = Ratio( "#gamma_{tight}/#gamma_{loose}", tight, loose )
		r.draw(0,2)
		c.SaveAs("plots/findFoId_%s_twoDistributions_%s.pdf"%(getSaveNameFromDatasets(filenames),''.join([s for s in fcut if s.isdigit()])) )

	return tight.Chi2Test( loose, option )


def optimiseFoId( filenames ):
	commonCut = "eventNumber%2==0"

	# default parameters
	v = {}
	v["chMax"] = 2.6
	v["nMax"] = 3.5
	v["nMaxRel"] = 0.04
	v["pMax"] = 1.3
	v["pMaxRel"] = 0.005

	v["chMin"] = -1
	v["nMin"] = -1
	v["pMin"] = -1

	latexText = {}
	latexText["chMax"] = "I_{#pm}^{max}"
	latexText["nMax"] = "I_{0}^{max}"
	latexText["nMaxRel"] = "I_{0,rel}^{max}"
	latexText["pMax"] = "I_{#gamma}^{max}"
	latexText["pMaxRel"] = "I_{#gamma,rel}^{max}"

	latexText["chMin"] = "I_{#pm}^{min}"
	latexText["nMin"] = "I_{0}^{min}"
	latexText["pMin"] = "I_{#gamma}^{min}"

	for abbr, text in latexText.iteritems():
		latexText[ abbr ] = text + " [GeV]"

	# tight histo does not change, so compute before
	tight = getHistoFromFiles( "met", "photonTree", filenames, commonCut )

	gr = ROOT.TGraph()
	gr.SetMaximum(1)
	gr.SetMinimum(0)

	vary = "chMax"

	for i, variable in enumerate(stepRange(2.6, 3.5, 20)):
		v[ vary ] = variable
		fcut = commonCut
		fcut += "&& photons[0].chargedIso<%s"%v["chMax"]
		fcut += "&& photons[0].neutralIso<%s+%s*photons[0].pt"%( v["nMax"], v["nMaxRel"] )
		fcut += "&& photons[0].photonIso <%s+%s*photons[0].pt"%( v["pMax"], v["pMaxRel"] )
		fcut += " &&( photons[0].chargedIso>%s && photons[0].neutralIso>%s && photons[0].photonIso>%s )"%(v["chMin"], v["nMin"], v["pMin"])
		chi2 = compareHists( filenames, tight, fcut, True )
		gr.SetPoint(i, variable, chi2 )


	foIdText = ROOT.TPaveText(.7,.85,.95,.95, "ndc")
	foIdText.SetFillColor(0)
	foIdText.SetBorderSize(0)
	vForTextInfo = v
	vForTextInfo[vary] = latexText[vary].replace(" [GeV]", "")
	foIdText.AddText( "I_{#pm} < %s"%vForTextInfo["chMax"] )
	foIdText.AddText( "I_{0} < %s + %s p_{T}"%(vForTextInfo["nMax"],vForTextInfo["nMaxRel"]) )
	foIdText.AddText( "I_{#gamma} < %s + %s p_{T}"%(vForTextInfo["pMax"],vForTextInfo["pMaxRel"]) )

	ROOT.gStyle.SetOptLogy(0)
	gr.SetTitle(";%s;p-Value"%latexText[vary] )

	c = ROOT.TCanvas()
	c.cd()
	gr.Draw("ap")
	foIdText.Draw()
	c.SaveAs( "plots/findFoId_%s_%s.pdf"%(getSaveNameFromDatasets(filenames), vary))
	ROOT.SetOwnership( c, False )
	del c








if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	opts = arguments.parse_args()

	optimiseFoId( opts.filenames )

