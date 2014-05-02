#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

info = PlotCaption()
info = ROOT.TLatex(0,.97, "#text{CMS Private Work - Wino-like }#tilde{#chi}^{0}_{1}#hspace{3cm}#, #sqrt{s}=#SI{8}{TeV}#, #geq1#ggamma,#geq2#text{jets}" )
info.SetNDC()
info.SetTextSize(1./15.8447)

ROOT.gStyle.SetCanvasDefH(1000)
ROOT.gStyle.SetCanvasDefW(2000)
ROOT.gStyle.SetPaperSize(14.65,50)

# Margins:
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.10)
ROOT.gStyle.SetPadRightMargin(0.02)

ROOT.gStyle.SetErrorX(0)

def drawSignalContamination( filename, xy, split ):
	import array
	label, unit, binning = readAxisConf("met")
	metBinning = array.array( "d", binning )

	gHisto = readHisto( filename, "gMet%s_%s"%xy ).Rebin(len(metBinning)-1, randomName(), metBinning )
	eHisto = readHisto( filename, "eMet%s_%s"%xy ).Rebin(len(metBinning)-1, randomName(), metBinning )
	fHisto = readHisto( filename, "fMet%s_%s"%xy ).Rebin(len(metBinning)-1, randomName(), metBinning )

	eHisto = applyFakeRateEWK( eHisto )
	eContamination = divideHistos( eHisto, gHisto )
	eContamination.SetLineColor(2)

	fContamination = divideHistos( fHisto, gHisto )
	fContamination.SetLineColor(4)

	totalContamination = divideHistos( addHistos( [eHisto, fHisto] ), gHisto )
	totalContamination.SetLineColor(1)
	totalContamination.SetMarkerColor(1)

	for h in [totalContamination, fContamination, eContamination]:
		h.Scale(100.) # in %
		h.GetYaxis().SetTitleOffset( 0.9 )
		h.SetTitle(";%s [#text{%s}];Signal Contamination [%s]"%("#met",unit, "%"))
		h.SetLineWidth(2)
		h.SetMarkerSize(0)
		h.SetLabelSize(1./15.8447, "xy")
		h.SetTitleSize(1./15.8447, "xy")

	mh = Multihisto()
	mh.setMinimum(0)
	if split:
		mh.addHisto( totalContamination, "Total", draw="pe" )
		mh.addHisto( fContamination, "b_{#text{signal}}^{#text{QCD}}/s", draw="hist" )
		mh.addHisto( eContamination, "b_{#text{signal}}^{e#rightarrow#gamma}/s", draw="hist" )
	else:
		mh.addHisto( totalContamination, "", draw="e0" )

	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(False)
	mh.Draw()
	totalContamination.Draw("same")
	info.Draw()
	saveName = "signalContamination_%s_%s_%s_%s"%(xy+(split,filename[0:-4] ))
	SaveAs(can, saveName )
	ROOT.gPad.SaveAs("/home/knut/master/documents/thesis/plots/%stex"%saveName )
	correctTiksPlot( "/home/knut/master/documents/thesis/plots/%stex"%saveName )





if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--xy", nargs=2, default = [1000,1220], type=int )
	arguments.add_argument("--split", action="store_true" )
	opts = arguments.parse_args()

	for inName in opts.filenames:
		drawSignalContamination( inName, tuple(opts.xy), opts.split )

