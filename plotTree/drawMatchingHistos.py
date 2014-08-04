#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *


def getPlots( filename ):
	plots = []
	# extract plot names only for the first file. Assume the files are similar
	f = ROOT.TFile( filename )
	for item in f.GetListOfKeys():
		name = item.GetName()
		if "match" in name:
			plots.append( ( name, readHisto( filename, name ) ) )
	return plots

def sumUpPlots( filenames, norm ):

	hists = getPlots( filenames[0] )

	for filename in filenames[1:]:
		for i in range(len(hists)):
			thisHists = getPlots( filename )
			hists[i][1].Add( thisHists[i][1] )

	datasetAffix = getSaveNameFromDatasets( filenames )

	Styles.tdrStyle2D()
	#ROOT.gStyle.SetNumberContours( 999 )

	#st = ROOT.gROOT.GetStyle("tdrStyle")
	#paperWidth = 14.65 #cm
	#st.SetPaperSize(paperWidth/2,50.)


	ROOT.gStyle.SetPadTopMargin(0.08)
	#st.SetPadRightMargin(0.1)
	for name, histo in hists:
		# deltaR
		histo.GetXaxis().SetLabelSize(0.06)
		histo.GetXaxis().SetTitleSize(0.06)
		histo.GetXaxis().SetTitleOffset(0.8)
		histo.GetXaxis().SetNdivisions(5,5,0,False)

		# rel pt
		histo.GetYaxis().SetLabelSize( histo.GetXaxis().GetLabelSize() )
		histo.GetYaxis().SetTitleSize( histo.GetXaxis().GetTitleSize() )
		histo.GetZaxis().SetLabelSize( histo.GetXaxis().GetLabelSize() )
		histo.GetZaxis().SetTitleSize( histo.GetXaxis().GetTitleSize() )


		if norm and histo.Integral():
			histo.Scale( 1./histo.Integral() )
		histo.Draw("colz")
		if "matchPhoton" in name and not name.endswith( "Pt" ):
			#histo.GetXaxis().SetTitle("$#Delta{R}$")
			#histo.GetYaxis().SetTitle("$p_{T, #text{jet}} / p_{T,#gamma}$")
			drCut = .2
			relPtCut = 3
			relPtMinCut = 0.8

			drCutLine = ROOT.TLine( drCut, relPtMinCut, drCut, relPtCut )
			relPtCutLine = ROOT.TLine( 0, relPtCut, drCut, relPtCut )
			relPtCutMinLine = ROOT.TLine( 0, relPtMinCut, drCut, relPtMinCut )

			for line in drCutLine, relPtCutLine, relPtCutMinLine:
				line.SetLineWidth(2)
				line.SetLineColor(2)
				line.Draw()

		if "matchGenPhoton" in name:
			histo.GetXaxis().SetTitle("#Delta{R}")
			histo.GetYaxis().SetTitle("(p_{T, #text{gen}}-p_T) / p_{T,#text{gen}}")

		pc1 = ROOT.TLatex(0,.95, "CMS Preliminary - 19.7fb^{-1}(8TeV)")
		if name == "matchPhotonToJet":
			pc1 = ROOT.TLatex(0,.95, "CMS Preliminary - 19.7fb^{-1}(8TeV) - #gamma_{tight}")
		if name == "matchPhotonJetToJet":
			pc1 = ROOT.TLatex(0,.95, "CMS Preliminary - 19.7fb^{-1}(8TeV) - #gamma_{loose}")
		for pc in [pc1]:
			pc.SetNDC()
			pc.SetTextSize(histo.GetXaxis().GetTitleSize()*.80)
			pc.SetTextFont(histo.GetXaxis().GetTitleFont())
			pc.Draw()


		SaveAs( ROOT.gPad, "%s_%s"%(name, datasetAffix) )
		#if name in ["matchPhotonToJet","matchPhotonJetToJet", "matchGenPhoton" ]:
		#	saveName = "/home/knut/master/documents/thesis/plots/%s_%s.tex"%(manipulateSaveName(name),shortName( filenames ))
		#	ROOT.gPad.SaveAs( saveName )
		#	correctTiksPlot( saveName )


if __name__ == "__main__":
	arguments = argparse.ArgumentParser( description="Calculate weighting "
			+"factors for QCD background estimation." )
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--sum", action="store_true" )
	arguments.add_argument("--norm", action="store_true" )
	opts = arguments.parse_args()


	if opts.sum:
		sumUpPlots( opts.filenames, opts.norm )
	else:
		for fileName in opts.filenames:
			sumUpPlots( [fileName], opts.norm )

