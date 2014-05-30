#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from predictions import *

def drawWeightHisto( weight2D, saveName, writeWeightFile=False, control=True ):
	regionString = "control" if control else "signal"
	# Draw the histograms
	info = PlotCaption(control=control, signal=not control,treeName="")
	info = ROOT.TLatex(0,.96, "CMS Private Work  - 19.7fb^{-1} #sqrt{s}=8TeV #geq1#gamma,#geq2jets #slash{E}_{T}<100GeV" )
	info.SetNDC()

	# Display the weight errors as 2D histograms.
	weight2D.SetName("weight2D")
	weightErrors = weight2D.Clone( "weightError" )
	weightRelErrors = weight2D.Clone( "weightRelError" )
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if weight2D.GetBinContent( i, j ):
				weightRelErrors.SetBinContent( i, j, weight2D.GetBinError( i, j )/weight2D.GetBinContent( i, j ) )

			weightErrors.SetBinContent( i, j, weight2D.GetBinError( i, j ) )

	# Draw histograms
	Styles.tdrStyle2D()
	ROOT.gStyle.SetPaintTextFormat("1.1f");

	can2D = ROOT.TCanvas()
	can2D.cd()
	can2D.SetLogz(1)

	for hist in weight2D, weightErrors, weightRelErrors:

		name = hist.GetName()
		if name == "weight2D":
			hist.GetZaxis().SetRangeUser(0.8,20)
			hist.GetZaxis().SetTitle("\qcdRatio")
		if name == "weightRelError":
			#hist.GetZaxis().SetRangeUser(0,1.2)
			hist.GetZaxis().SetTitle("\qcdRatioError")

		hist.Draw("colz")
		info.Draw()
		SaveAs(can2D, "qcd_preWeight_%s_%s_%s"%(saveName,name,regionString) )

	Styles.tdrStyle()

	if writeWeightFile:
		weightFile = ROOT.TFile( "qcdWeight.root", "recreate" )
		weightFile.cd()
		weight2D.SetName("qcdWeight")
		weight2D.Write(ROOT.TObject.kOverwrite)
		weightFile.Close()



def drawTwoHists( gHist, fHist, sHist, saveName, minmax ):
	for bin in range( sHist.GetNbinsX()+2 ):
		sHist.SetBinError( bin, sHist.GetBinContent(bin) )
		sHist.SetBinContent( bin, fHist.GetBinContent(bin) )

	sHist.SetFillColor( sHist.GetLineColor() )
	sHist.SetLineColor( sHist.GetLineColor() )
	sHist.SetFillStyle(3254)
	sHist.SetMarkerSize(0)
	fHist.SetLineColor(2)
	sHist.SetLineColor(2)


	for h in gHist, fHist, sHist:
		from inheritRoot import H1F
		h.__class__ = H1F
		h.MergeOverflow()
		h.Scale(1., "width")

	muhisto = Multihisto()
	muhisto.addHisto( gHist, "Simulation", draw="hist e" )
	muhisto.addHisto( fHist, "Prediction", draw="hist")
	muhisto.addHisto( sHist, "#sigma_{w}", draw="e2")

	can = ROOT.TCanvas(randomName(), "", 1000, 1200)
	can.cd()
	muhisto.Draw()
	text = ROOT.TLatex(.1,.965, "%i #leq p_{T^{*}} < %i, %i #leq H_{T} < %i"%minmax)
	text.SetNDC()
	text.Draw()

	from myRatio import Ratio
	r = Ratio( "Sim./Pred.", gHist, fHist )
	r.draw(0,2)

	can.SaveAs(saveName+".pdf")

	ROOT.SetOwnership( can, False )
	del can



def drawSignleClosure( filenames, plot, commonCut, info ):
	from ROOT import TH3F
	import array
	ptBins = readAxisConf( "photons[0].ptJet()" )[2]
	htBins = readAxisConf( "ht" )[2]
	metBins = readAxisConf( "met" )[2]

	ptBinsAr = array.array( 'd', ptBins )
	htBinsAr = array.array( 'd', htBins )
	metBinsAr = array.array( 'd', metBins )

	gname = randomName()
	fname = randomName()
	sname = randomName()

	goHist3d = TH3F( gname, "3d prediction", len(ptBins)-1, ptBinsAr, len(htBins)-1, htBinsAr, len(metBins)-1, metBinsAr )
	foHist3d = TH3F( fname, "3d prediction", len(ptBins)-1, ptBinsAr, len(htBins)-1, htBinsAr, len(metBins)-1, metBinsAr )
	sHist3d = TH3F( sname, "3d prediction", len(ptBins)-1, ptBinsAr, len(htBins)-1, htBinsAr, len(metBins)-1, metBinsAr )

	for h in goHist3d, foHist3d, sHist3d:
		h.Sumw2()

	for filename in filenames:
		goTree = readTree( filename, "photonTree" )
		goTree.Draw( "met:ht:photons[0].ptJet()>>+%s"%gname, "weight", "goff" )

		foTree = readTree( filename, "photonJetTree" )
		foTree.AddFriend( "foWeights", filename )
		foTree.Draw( "met:ht:photons[0].ptJet()>>+%s"%fname, "weight*w_qcd", "goff" )
		foTree.Draw( "met:ht:photons[0].ptJet()>>+%s"%sname, "weight*w_qcd_error", "goff" )

	for ptIndex in range(0, len(ptBins)-1):
		for htIndex in range(0, len(htBins)-1):
			gHist = goHist3d.ProjectionZ( randomName(), ptIndex+1, ptIndex+1, htIndex+1,htIndex+1, "e")
			fHist = foHist3d.ProjectionZ( randomName(), ptIndex+1, ptIndex+1, htIndex+1,htIndex+1, "e")
			sHist = sHist3d.ProjectionZ( randomName(), ptIndex+1, ptIndex+1, htIndex+1,htIndex+1, "e")
			drawTwoHists( gHist, fHist, sHist, "test_%02d_%02d"%(ptIndex,htIndex), (ptBins[ptIndex], ptBins[ptIndex+1], htBins[htIndex], htBins[htIndex+1]) )


def drawClosure( filenames, predFilenames, plot, commonCut, infoText, additionalLabel="", modifyEmptyBins=True ):
	infoText = ROOT.TLatex(0,.96, "CMS Private Work 19.7fb^{-1} #sqrt{s}=8TeV #geq1#gamma,#geq2jets" )
	infoText.SetNDC()
	infoText.SetTextSize(0.045)

	gHist = getHists( filenames, plot, commonCut )
	fHist, sysHist = predictionHistos( predFilenames, plot, commonCut, modifyEmptyBins )
	sysHist.SetFillStyle(3254)
	sysHist.SetFillColor( sysHist.GetLineColor() )

	for bin in range(gHist.FindBin(101), gHist.GetNbinsX()+1):
		w = gHist.GetBinWidth(bin)
		print gHist.GetBinLowEdge(bin)
		print gHist.GetBinContent(bin)*w, "±", gHist.GetBinError(bin)*w
		print fHist.GetBinContent(bin)*w, "±", fHist.GetBinError(bin)*w, "±", sysHist.GetBinContent(bin)*w

	signalAbbrs = mergeDatasetAbbr( [ getDatasetAbbr(x) for x in filenames ] )

	muhisto = Multihisto()
	muhisto.leg.SetHeader( ",".join( [ datasetToLatex(x) for x in signalAbbrs ] ) )
	muhisto.addHisto( gHist, "Simulation", draw="hist e0" )
	muhisto.addHisto( fHist, "Prediction", draw="hist e")
	muhisto.addHisto( sysHist, "#sigma_{w}", draw="e2")
	muhisto.addHisto( gHist, "", draw="hist e" ) # add a second time to draw on top

	can = ROOT.TCanvas("", "", 1000, 1200)
	can.cd()

	muhisto.Draw()

	from myRatio import Ratio
	r = Ratio( "Sim./Pred.", gHist, fHist, sysHist )
	r.draw(0,2.4)
	infoText.Draw()
	SaveAs(can, "qcdClosure_%s_%s"%("".join(signalAbbrs)+additionalLabel, plot) )

	# Since root is too stupid to clear the canvas before python is ending, clean
	# the canvas yourself
	ROOT.SetOwnership( can, False )
	del can


def qcdClosure( filenames, plots ):
	signalCut = "met>=100"
	controlCut = "!(%s)"%signalCut
	commonCut = "!@electrons.size() && !@muons.size()"

	# Definition of labels
	infoControl = PlotCaption( control=True )
	infoSignal = PlotCaption( signal=True )
	info = PlotCaption()

	weights = getMixedWeigthHisto( filenames, filenames, commonCut )
	attachWeightsToFiles( filenames, weights, "foWeights" )
	#drawWeightHisto( weights, getSaveNameFromDatasets(filenames) )

	for plot in plots:
		if plot != "met":
			drawClosure( filenames, fenames, plot, commonCut+"&&"+signalCut, infoSignal, "_signal" )
			drawClosure( filenames, filenames, plot, commonCut+"&&"+controlCut, infoControl, "_control" )
		drawClosure( filenames, filenames, plot, commonCut, info, modifyEmptyBins=True )
		#drawSignleClosure( filenames, plot, commonCut, info )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("-p", "--plots", nargs="+", default=["met"] )
	opts = arguments.parse_args()

	if opts.plots == ["all"]:
		opts.plots = [ "met", "ht", "photons[0].ptJet()","nGoodJets" ]

	qcdClosure( opts.filenames, opts.plots )

