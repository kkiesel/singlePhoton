#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
from sys import stdout

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
			stdout.write( "\r%s: %s / %s"%(fileName, event.GetReadEntry(), event.GetEntries() ) )
			stdout.flush()

		b = h_weight.FindBin( event.photons.at(0).ptJet(), event.ht )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()
	print

	f = ROOT.TFile( fileName, "update" )
	f.cd()
	weightTree.Write()
	f.Close()

def getMixedWeigthHisto( filenames, predFilenames, commonCut, control=True ):
	"""Calculate #photons/#photonFakes in bins of photons.ptJet and a second
	(global) variable.

	filenames: files containing photons
	predFilenames: files containing fakes
	"""

	regionCut = "met<100" if control else "met>=100"

	xVar = "photons[0].ptJet()"
	yVar = "ht"
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

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

	# Set the weight and error for empty bins above the diagonal to one.
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if not weight2D.GetBinContent( i, j ) and weight2D.GetXaxis().GetBinLowEdge(i) < weight2D.GetYaxis().GetBinUpEdge(j):
				weight2D.SetBinContent( i, j, 0 )
				weight2D.SetBinError( i, j, 0 )

	drawWeightHisto( weight2D, numerator, denominator, control, shortName( filenames ) )
	return weight2D

def drawWeightHisto( weight2D, numerator, denominator, control, saveName ):
	regionString = "control" if control else "signal"
	# Draw the histograms
	info = PlotCaption(control=control, signal=not control,treeName="")
	info = ROOT.TLatex(0,.96, "#text{CMS Private Work  -  }#SI{19.8}{fb^{-1}}#, #sqrt{s}=#SI{8}{TeV}#, #geq1#gamma,#geq2#text{jets}#,, #met<#SI{100}{GeV}" )
	info.SetNDC()
	info.SetTextSize(0.07930341347505648/1.63823)


	for h in weight2D, numerator, denominator:
		h.SetTitle(";#pt#text{ [GeV]};H_{T}#text{ [GeV]}")

	# Display the weight errors as 2D histograms.
	weight2D1 = divideHistos( numerator, denominator )
	weightErrors = weight2D.Clone( randomName() )
	weightErrors1 = weight2D.Clone( randomName() )
	weightRelErrors = weight2D.Clone( randomName() )
	weightRelErrors1 = weight2D.Clone( randomName() )
	for i in range( weight2D.GetXaxis().GetNbins()+1 ):
		for j in range( weight2D.GetYaxis().GetNbins()+1 ):
			if weight2D1.GetBinContent( i, j ):
				weightRelErrors.SetBinContent( i, j, weight2D.GetBinError( i, j )/weight2D.GetBinContent( i, j ) )
				weightRelErrors1.SetBinContent( i, j, weight2D1.GetBinError( i, j )/weight2D1.GetBinContent( i, j ) )

			weightErrors.SetBinContent( i, j, weight2D.GetBinError( i, j ) )
			weightErrors1.SetBinContent( i, j, weight2D1.GetBinError( i, j ) )

	# Draw histograms
	Styles.tdrStyle2D()
	ROOT.gStyle.SetPaintTextFormat("1.1f");
	ROOT.gStyle.SetPaperSize(12,50.)
	ROOT.gStyle.SetPadTopMargin(0.06)

	can2D = ROOT.TCanvas()
	can2D.cd()
	can2D.SetLogz(0)
	for hist, name in [
			(numerator,"numerator"),
			(denominator, "denominator")]:
		#hist.Draw("colz text")
		hist.Draw("colz")
		info.Draw()
		SaveAs(can2D, "qcd_preWeight_%s_%s_%s"%(saveName, name,regionString) )

	for hist, hist1, name in [
			(weight2D,weight2D1, "weight2D"),
			(weightErrors,weightErrors1, "weightError"),
			(weightRelErrors,weightRelErrors1, "weightRelError") ]:
		if name == "weight2D":
			hist1.GetZaxis().SetRangeUser(0,5)
			hist1.GetZaxis().SetTitle("\qcdRatio")
		if name == "weightRelError":
			hist1.GetZaxis().SetRangeUser(0,1.2)
			hist1.GetZaxis().SetTitle("\qcdRatioError")

		if name == "weight2D" or name == "weightRelError":
			for ax in hist1.GetXaxis(), hist1.GetYaxis(), hist1.GetZaxis():
				ax.SetLabelSize(1./20.6577)
				ax.SetTitleSize(1./20.6577)
				ax.SetTitleOffset(1.2)


		hist1.Draw("colz")
		#hist.Draw("same text")

		# Draw diagonal:
		if True:
			x0 = hist1.GetXaxis().GetBinLowEdge(1)
			x1 = hist1.GetXaxis().GetBinLowEdge(hist1.GetNbinsX()+1)
			y0 = hist1.GetYaxis().GetBinLowEdge(1)
			y1 = hist1.GetYaxis().GetBinLowEdge(hist1.GetNbinsY()+1)
			# create diagonal
			l = ROOT.TLine( max(x0,y0), max(x0,y0), min(x1,y1), min(x1,y1) )
			#l.Draw()

			# a second function
			# function f
			f = lambda x: 2*x-150
			# inverse of this function
			F = lambda x: 0.5*(x+150)
			points = [ x0, f(x0) ] if x0 > y0 else [ F(y0), y0 ]
			points += [ x1, f(x1) ] if x1 > y1 else [ F(y1), y1 ]
			l2 = ROOT.TLine( *points )
			#l2.Draw()

		info.Draw()
		SaveAs(can2D, "qcd_preWeight_%s_%s_%s"%(saveName,name,regionString) )
		if name == "weight2D" or name == "weightRelError":
			can2D.SaveAs("~/master/documents/thesis/plots/qcd_%s_%s_%s.tex"%(saveName,name,regionString) )

	Styles.tdrStyle()

	if control:
		weightFile = ROOT.TFile( "qcdWeight.root", "recreate" )
		weightFile.cd()
		weight2D.SetName("qcdWeight")
		weight2D.Write()
		weightFile.Close()

def photonHisto( filenames, plot, cut, modifyEmptyBins ):
	gHist = None
	for filename in filenames:
		gTree = readTree( filename, "photonTree" )
		hist = getHisto( gTree, plot, cut=cut, color=1, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		if gHist:
			gHist.Add( hist )
		else:
			gHist = hist
	return gHist

def predictionHistos( filenames, plot, cut, modifyEmptyBins ):
	fHist, sysHist = None, None

	for filename in filenames:
		fTree = readTree( filename, "photonJetTree" )
		fTree.AddFriend( "foWeights", filename )

		# Prediction + statistical uncertainty coming from the loose control region
		hist = getHisto( fTree, plot, weight="weight*w_qcd", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )

		# This histogram contains the uncertainty due to the weight for small met
		# Here, the bin content is the uncertainty, the bin error is meaningless.
		sHist = getHisto( fTree, plot, weight="weight*w_qcd_error", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )

		# You can also do an approximation to calculate the error of the weight:
		doWeightApproximation = False
		if doWeightApproximation:
			sHistUp = getHisto( fTree, plot, weight="weight*(w_qcd_error+w_qcd)",
					cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
			sHistDown = getHisto( fTree, plot, weight="weight*(w_qcd-w_qcd_error)",
					cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
			for bin in range( sHistUp.GetNbinsX()+2 ):
				sHistUp.SetBinContent( bin, 0.5*(sHistUp.GetBinContent(bin)-sHistDown.GetBinContent(bin)))
			sHist = sHistUp
		# Approximation end ###################################################

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

	return fHist, sysHist


def drawClosure( filenames, predFilenames, plot, commonCut, infoText, additionalLabel, modifyEmptyBins=True ):
	infoText = ROOT.TLatex(0,.96, "#text{CMS Private Work#hspace{5cm}   }#SI{19.8}{fb^{-1}}#, #sqrt{s}=#SI{8}{TeV}#, #geq1#gamma,#geq2#text{jets}" )
	infoText.SetNDC()
	infoText.SetTextSize(0.07930341347505648/1.63823/1.52236)

	gHist = photonHisto( filenames, plot, commonCut, modifyEmptyBins )
	fHist, sysHist = predictionHistos( predFilenames, plot, commonCut, modifyEmptyBins )

	for bin in range( sysHist.GetNbinsX()+2 ):
		sysHist.SetBinError( bin, sysHist.GetBinContent(bin) )
		sysHist.SetBinContent( bin, fHist.GetBinContent(bin) )


	signalAbbrs = mergeDatasetAbbr( [ getDatasetAbbr(x) for x in filenames ] )
	predAbbrs =  mergeDatasetAbbr( [ getDatasetAbbr(x) for x in predFilenames ] )

	for h in gHist, fHist,sysHist:
		h.GetXaxis().SetTitle("#met#text{ [GeV]}")
		for a in h.GetXaxis(), h.GetYaxis():
			a.SetLabelSize(1./31.4485)
			a.SetTitleSize(1./31.4485)
		h.GetXaxis().SetTitleOffset(1.1)
		h.GetYaxis().SetTitleOffset(1.4)

	muhisto = Multihisto()
	muhisto.leg.SetX1(.63)
	muhisto.leg.SetX2(.95)
	#muhisto.leg.SetHeader( ",".join( [ datasetToLatex(x) for x in signalAbbrs ] ) )
	muhisto.addHisto( gHist, "Simulation", draw="hist e" )
	muhisto.addHisto( fHist, "Prediction", draw="hist")
	muhisto.addHisto( sysHist, "#sigma_{#met<#SI{100}{GeV}}", draw="e2")
	muhisto.addHisto( gHist, "", draw="hist e" ) # add a second time to draw on top

	ROOT.gStyle.SetPaperSize(14.6,50.)
	ROOT.gStyle.SetPadTopMargin(0.05)
	ROOT.gStyle.SetPadRightMargin(0.02)
	ROOT.gStyle.SetPadLeftMargin(0.09)

	can = ROOT.TCanvas()
	can.cd()

	muhisto.Draw()

	fHistSumError = fHist.Clone( randomName() )
	for bin in range( fHistSumError.GetNbinsX()+2 ):
		fHistSumError.SetBinError( bin, sqrt( fHistSumError.GetBinError(bin)**2+sysHist.GetBinError(bin)**2 ) )
		gHist.SetBinError( bin, sqrt( gHist.GetBinError(bin)**2 + fHist.GetBinError(bin)**2 ) )

	from myRatio import Ratio
	#r = Ratio( "#gamma/#gamma_{pred}", gHist, fHistSumError )
	r = Ratio( "Sim./Pred.", gHist, sysHist )
	r.draw(0,2)
	infoText.Draw()
	SaveAs(can, "qcdClosure_%s_%s"%("".join(signalAbbrs)+additionalLabel, plot) )
	can.SaveAs("~/master/documents/thesis/plots/qcdClosure_%s_%s.tex"%(''.join(signalAbbrs), plot) )
	correctTiksPlot( "/home/knut/master/documents/thesis/plots/qcdClosure_%s_%s.tex"%(''.join(signalAbbrs), plot) )


	# Since root is too stupid to clear the canvas before python is ending, clean
	# the canvas yourself
	ROOT.SetOwnership( can, False )
	del can



def qcdClosure( filenames, predFilenames, plots ):
	signalCut = "met>=100"
	controlCut = "!(%s)"%signalCut
	commonCut = "1"
	#commonCut += " && photons[0].isGen(2)"
	#commonCut = "!photons[0].isGen(0) && !photons[0].isGen(1)"

	# Definition of labels
	infoControl = PlotCaption( control=True )
	infoSignal = PlotCaption( signal=True )
	info = PlotCaption()

	weights = getMixedWeigthHisto( filenames, predFilenames, commonCut )
	getMixedWeigthHisto( filenames, predFilenames, commonCut, False )

	for fileName in predFilenames:
		foTree = readTree( fileName, "photonJetTree" )
		writeWeight2DToFile( fileName, foTree, weights, "foWeights" )

	for plot in plots:
		if plot != "met":
			drawClosure( filenames, predFilenames, plot, commonCut+"&&"+signalCut, infoSignal, "_signal" )
			drawClosure( filenames, predFilenames, plot, commonCut+"&&"+controlCut, infoControl, "_control" )
		drawClosure( filenames, predFilenames, plot, commonCut, info, "", True )

if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("-p", "--prediction", nargs="+", default=[], type=isValidFile )
	arguments.add_argument("--plot", nargs="+", default=["met"] )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "photons[0].ptJet()","nGoodJets" ]
	if not opts.prediction:
		opts.prediction = opts.filenames

	qcdClosure( opts.filenames, opts.prediction, opts.plot )

