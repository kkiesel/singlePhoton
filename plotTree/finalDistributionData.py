#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

# definition of an Infix operator class
# this recipe also works in jython
# calling sequence for the infix is either:
#  x |op| y
# or:
# x <<op>> y

class Infix:
    def __init__(self, function):
        self.function = function
    def __ror__(self, other):
        return Infix(lambda x, self=self, other=other: self.function(other, x))
    def __or__(self, other):
        return self.function(other)
    def __rlshift__(self, other):
        return Infix(lambda x, self=self, other=other: self.function(other, x))
    def __rshift__(self, other):
        return self.function(other)
    def __call__(self, value1, value2):
        return self.function(value1, value2)

# quadratic addition, usage: 3 |qPlus| 4
qPlus = Infix( lambda x,y: sqrt(x**2+y**2) )




def printDeviations( h1, h2 ):
	print "Deviations:"
	for bin in range( h1.GetNbinsX()+2):
		deviation = h1.GetBinContent(bin) - h2.GetBinContent(bin)
		error = sqrt( h1.GetBinError(bin)**2 + h2.GetBinError(bin)**2)
		if error:
			print h1.GetBinLowEdge(bin), deviation / error


def getMixedWeigthHisto( filenames, predFilenames, commonCut, control=True, fillEmptyBins=False ):
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
		den = createHistoFromTree2D( foTree, yVar+":"+xVar, "weight*( %s && %s && !photons.isGen(2) && @photons.size() )"%(regionCut, commonCut), xbinning, ybinning )
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
	weightTree.Branch( "w_qcd", weight, "w_qcd/D" )
	weightTree.Branch( "w_qcd_error", weight_error, "w_qcd_error/D" )

	from sys import stdout
	for event in tree:
		if not event.GetReadEntry()%10000:
			stdout.write( "\r%s / %s"%(event.GetReadEntry(), event.GetEntries() ) )
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

def drawWeightHisto( weight, control=True ):
	regionString = "control" if control else "signal"
	# Draw the histograms
	info = PlotCaption(control=control, signal=not control,treeName="")
	info = ROOT.TLatex(0,.96, "#text{CMS Private Work  #hspace{2cm}  }#SI{19.8}{fb^{-1}}#, #sqrt{s}=#SI{8}{TeV}#, #geq1#gamma,#geq2#text{jets}#,, #met<#SI{100}{GeV}" )
	info.SetNDC()
	info.SetTextSize(0.07930341347505648/1.63823/1.16666)


	weight.SetTitle(";#pt#text{ [GeV]};H_{T}#text{ [GeV]}")

	# Display the weight errors as 2D histograms.
	weightErrors = weight.Clone( randomName() )
	weightRelErrors = weight.Clone( randomName() )
	for i in range( weight.GetXaxis().GetNbins()+1 ):
		for j in range( weight.GetYaxis().GetNbins()+1 ):

			weightErrors.SetBinContent( i, j, weight.GetBinError( i, j ) )
			if weight.GetBinContent( i, j ):
				weightRelErrors.SetBinContent( i, j, weight.GetBinError( i, j )/weight.GetBinContent( i, j ) )

	# Draw histograms
	Styles.tdrStyle2D()
	ROOT.gStyle.SetPaintTextFormat("1.1f");
	ROOT.gStyle.SetPaperSize(14,50.)
	ROOT.gStyle.SetPadTopMargin(0.06)

	can2D = ROOT.TCanvas()
	can2D.cd()
	can2D.SetLogz(0)

	for hist, zTitle, name in [
			(weight, "#ggamma / #fgamma", "weight"),
#			(weightErrors,"weighterr", "weightError"),
			(weightRelErrors,"#frac{#sigma_{#ggamma/#fgamma}}{#ggamma/#fgamma}", "weightRelError") ]:
		if name == "weight":
			hist.GetZaxis().SetRangeUser(0,5)
		if name == "weightRelError":
			hist.GetZaxis().SetRangeUser(0,1.2)
		hist.GetZaxis().SetTitle( zTitle )

		for ax in hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis():
			ax.SetLabelSize(1./20.6577/1.16667)
			ax.SetTitleSize(1./20.6577/1.16667)
		hist.GetXaxis().SetTitleOffset(1)
		hist.GetYaxis().SetTitleOffset(1.4)
		hist.GetZaxis().SetTitleOffset(1.1)

		hist.Draw("colz")

		info.Draw()
		SaveAs(can2D, "qcd_preWeight_data_%s_%s"%(name,regionString) )
		if name == "weight" or name == "weightRelError":
			can2D.SaveAs("~/master/documents/thesis/plots/qcd_data_%s_%s.tex"%(name,regionString) )

	Styles.tdrStyle()

def getHists( filenames, cut="1", plot="met", treeName="photonTree" ):
	endHist = None
	for filename in filenames:
		tree = readTree( filename, treeName )
		hist = getHisto( tree, plot, color=1, fillEmptyBins=not ("PhotonHad" in filename), cut=cut )

		if endHist:
			endHist.Add( hist )
		else:
			endHist = hist

	return endHist

def applyAddUncertainty( hist, uncert ):
	for bin in range( hist.GetNbinsX()+1 ):
		hist.SetBinError( bin, uncert * hist.GetBinContent(bin) )
	return hist


def qcdPredictionHistos( filenames, plot, cut, modifyEmptyBins ):
	fHist, sysHist, sysEmptyBin = None, None, None
	cut+=" && !photons.isGen(2) && @photons.size()"
	for filename in filenames:
		fTree = readTree( filename, "photonJetTree" )
		fTree.AddFriend( "foWeights", filename )

		hist = getHisto( fTree, plot, weight="weight*w_qcd", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		sHist = getHisto( fTree, plot, weight="weight*w_qcd_error", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		sHistEmptyBin = getHisto( fTree, plot, weight="weight", cut="w_qcd*(%s)+(w_qcd<0.0001)*(%s)"%(cut,cut), color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
		sHistEmptyBin.Add( hist, -1. ) # get difference to normal prediction

		if fHist:
			fHist.Add( hist )
			sysHist.Add( sHist )
			sysEmptyBin.Add( sHistEmptyBin )
		else:
			fHist = hist
			sysHist = sHist
			sysEmptyBin = sHistEmptyBin

	sysHist.SetFillColor( sysHist.GetLineColor() )
	sysHist.SetLineColor( sysHist.GetLineColor() )
	sysHist.SetFillStyle(3254)
	sysHist.SetMarkerSize(0)

	fHist.SetLineColor(7)

	return fHist, sysHist, sysEmptyBin

def printBin( histo, error=False, scale=1., roundInt=False ):
	# roundInt is necessary since ROOT has a low precission, leading to values
	# close to an integer (eg 23.000048)

	out = ""

	import array
	minimalMet=100
	hClone = histo.Clone( randomName() )
	hClone.Scale( scale )
	for bin in range(1, hClone.GetNbinsX()+1):
		if bin<hClone.FindBin(100):
			continue

		if error:
			out += '%s '%(hClone.GetBinError(bin)*hClone.GetBinWidth(bin))
		else:
			if roundInt:
				out += '%s '%int(round(hClone.GetBinContent(bin)*hClone.GetBinWidth(bin)))
			else:
				out += "%s "%(hClone.GetBinContent(bin)*hClone.GetBinWidth(bin))
	return out[0:-1]+"\n"

def getArrayOfBins( histo, error=False, scale=1., roundInt=False ):
	# roundInt is necessary since ROOT has a low precission, leading to values
	# close to an integer (eg 23.000048)

	import array
	minimalMet=100
	hClone = histo.Clone( randomName() )
	hClone.Scale( scale )

	out = []

	for bin in range(1, hClone.GetNbinsX()+1):
		if bin<hClone.FindBin(100):
			continue

		if error:
			out.append( hClone.GetBinError(bin)*hClone.GetBinWidth(bin) )
		else:
			if roundInt:
				out.append( int(round(hClone.GetBinContent(bin)*hClone.GetBinWidth(bin))) )
			else:
				out.append( hClone.GetBinContent(bin)*hClone.GetBinWidth(bin) )
	return out


def getIntAndError( h, bin ):
	err = ROOT.Double()
	con = h.IntegralAndError( bin, bin, err, "width" )
	return con, err


def writeFinalTableLong( versionData, fakeRate, fakeRateStatError, fakeRateSysError, fakeRateSysErrorOwn,
		isrUncertaintyZ, isrUncertaintyW, isrUncertaintyT,
		dataHist, fgammaHist, fgammaWeightError, fgammaEmptyBinError,
		egammaHist, fsrZ, fsrW, fsrT ):


	totalISR = addHistos( [fsrZ, fsrW, fsrT] )

	# Print out events bin by bin
	binList = [100.0,]
	for bin in range( dataHist.GetNbinsX()+1 ):
		xval = dataHist.GetBinLowEdge(bin)
		if xval > 100:
			binList.append( xval )

	import time
	dataCardString = """
# This file contains event yield and uncertainties for data and simulated event points.
# Produced by:
#    Knut Kiesel
# on:
#    %s
# with input version:
#    %s
# There will be a common section for all signal points, folllowed by the information
# for each signal point. The numbering of the signal points is arbitrary.
# For bin edges, the ROOT conversion is used (lower bin edge included, upper bin
# edge excluded).
"""%( time.strftime("%Y-%m-%d %H:%M:%S"), versionData )

	commonInformation = """
lumi = 3932 #pb
nMetBins = %s
""" % len(binList)

	binListInf = binList + ['infinity']

	for i in range(len(binList)):
		commonInformation += "bin %s = %s to %s\n"%(i,binListInf[i],binListInf[i+1])

	dataCardString += "\n##################################\n"
	dataCardString += commonInformation
	dataCardString += "\n##################################\n"

	data = getArrayOfBins( dataHist, roundInt=True )
	qcd = getArrayOfBins( fgammaHist )
	qcd_stat = getArrayOfBins( fgammaHist, True )
	qcd_weight = getArrayOfBins( fgammaWeightError )
	qcd_emptyBin = getArrayOfBins( fgammaEmptyBinError, roundInt=True )

	ewk = getArrayOfBins( egammaHist )
	ewk_stat = getArrayOfBins( egammaHist, True )
	ewk_f_stat = getArrayOfBins( egammaHist, scale=fakeRateStatError/fakeRate )
	ewk_f_syst = getArrayOfBins( egammaHist, scale=fakeRateSysError/fakeRate )
	ewk_transfer = getArrayOfBins( egammaHist, scale=fakeRateSysErrorOwn/fakeRate )

	isr = getArrayOfBins( totalISR )
	isr_stat = getArrayOfBins( totalISR, True )
	isr_sys = getArrayOfBins( totalISR, scale=.5 )

	isrW = getArrayOfBins( fsrW )
	isrZ = getArrayOfBins( fsrZ )
	isrT = getArrayOfBins( fsrT )

	isrW_stat = getArrayOfBins( fsrW, True )
	isrZ_stat = getArrayOfBins( fsrZ, True )
	isrT_stat = getArrayOfBins( fsrT, True )

	isrW_syst = getArrayOfBins( fsrW, scale=isrUncertaintyW )
	isrZ_syst = getArrayOfBins( fsrZ, scale=isrUncertaintyZ )
	isrT_syst = getArrayOfBins( fsrT, scale=isrUncertaintyT )

	for name, array in [
		( "selected        = %s", data ),
		( "QCD background  = %s", qcd ),
		( "QCD stat uncert = %s", qcd_stat ),
		( "QCD syst uncert = %s", qcd_weight ),
		( "QCD syst uncert = %s", qcd_emptyBin ),
		( "EWK background  = %s", ewk ),
		( "EWK stat uncert = %s", ewk_stat ),
		( "EWK stat uncert = %s", ewk_f_stat ),
		( "EWK syst uncert = %s", ewk_f_syst ),
		( "EWK syst uncert = %s", ewk_transfer ),
		( "ISR background  = %s", isr ),
		( "ISR stat uncert = %s", isr_stat ),
		( "ISR syst uncert = %s", isr_sys ) ]:
		dataCardString += name%( ' '.join(map(str, array ) ) ) + "\n"

	print dataCardString

	dataCardFile = open("eventYieldData-%s.txt"%time.strftime("%Y-%m-%d"), "w")
	dataCardFile.write( dataCardString )
	dataCardFile.close()

	rawTable = []
	rawTable.append([])
	rawTable[-1].append( "QCD" )
	rawTable[-1].append( "EWK" )
	rawTable[-1].append( "$\gamma$W" )
	rawTable[-1].append( "$t\\bar{t}$" )
	rawTable[-1].append( "$\gamma$Z" )
	rawTable[-1].append( "\\hline\nSum" )
	rawTable[-1].append( "Data" )

	binsToPrint = range(len(qcd) )
	binsToPrint = [0,5]

	for Bin in binsToPrint:
		rawTable.append([])
		rawTable[-1].append( "$%.2f \pm %.2f \pm %.2f \pm %i$"%( qcd[Bin], qcd_stat[Bin], qcd_weight[Bin], qcd_emptyBin[Bin] ) )
		rawTable[-1].append( "$%.2f \pm %.2f \pm %.2f \pm %.2f \pm %.2f$"%( ewk[Bin], ewk_stat[Bin], ewk_f_stat[Bin], ewk_f_syst[Bin], ewk_transfer[Bin] ) )
		rawTable[-1].append( "$%.2f \pm %.2f \pm %.2f$"%( isrW[Bin], isrW_stat[Bin], isrW_syst[Bin] ) )
		rawTable[-1].append( "$%.2f \pm %.2f \pm %.2f$"%( isrT[Bin], isrT_stat[Bin], isrT_syst[Bin] ) )
		rawTable[-1].append( "$%.2f \pm %.2f \pm %.2f$"%( isrZ[Bin], isrZ_stat[Bin], isrZ_syst[Bin] ) )
		rawTable[-1].append( "$%.2f \pm %.2f \pm %.2f$"%( qcd[Bin]+ewk[Bin]+isr[Bin], qcd_stat[Bin] |qPlus| ewk_stat[Bin] |qPlus| ewk_f_stat[Bin],
			qcd_weight[Bin]|qPlus|qcd_emptyBin[Bin]|qPlus| ewk_f_syst[Bin] |qPlus| ewk_transfer[Bin] |qPlus| isr_stat[Bin] |qPlus| isr_sys[Bin] ) )
		rawTable[-1].append( "$%s$"%data[Bin] )

	print "\\begin{tabular}{l|%s}"%(' '.join( ['l']*len(qcd) ) )
	print "\\hline\\hline"
	for Bin in binsToPrint:
		thisBin = 11+Bin
		print ' & $%i\leq\:\met[\si{GeV}] < %i$'%(dataHist.GetBinLowEdge(thisBin),dataHist.GetBinLowEdge(thisBin)+dataHist.GetBinWidth(thisBin) ),
	print " \\\\"
	print "\\hline"


	for line in zip(*rawTable):
		print ' & '.join(line) + " \\\\"
	print "\\hline\\hline"
	print "\\end{tabular}"



def finalDistributionData( plot ):

	# Some definitions:

	# Electroweak fake-rate
	fakeRate            = 1.48/100
	fakeRateStatError   = 0.05/100
	fakeRateSysError    = 0.08/100
	fakeRateSysErrorOwn = fakeRate * .5

	# additional ISR uncertainty
	isrUncertaintyZ = 0.7
	isrUncertaintyW = 0.5
	isrUncertaintyT = 0.5

	# Sample names
	wg1 = "slimWGamma_50_130_V02.%s_tree.root"%44
	wg2 = "slimWGamma_130_inf_V02.%s_tree.root"%44
	tg = "slimTTGamma_V02.%s_tree.root"%44
	zgn = "slimZGammaNuNu_V02.%s_tree.root"%43
	zgl = "slimZGammaLL_V02.%s_tree.root"%43

	versionData=44
	data = [ "PhotonHad%s_V02.%s_tree.root"%(x,versionData) for x in [ "A", "B", "C", "D" ] ]

	commonCut = "!@electrons.size() && !@muons.size()"

	blind = True
	additionalCut = "&& met < 100" if blind else ""

	# Compute the weights:
	#weight2D = getMixedWeigthHisto( data, data, commonCut, control=True, fillEmptyBins=False )
	#drawWeightHisto( weight2D )
	#for filename in data:
	#	fTree = readTree( filename, "photonJetTree" )
	#	writeWeight2DToFile( filename, fTree, weight2D, "foWeights" )

	# Get Histograms
	dataHist = getHists( data, commonCut+additionalCut )

	fgammaHist, fgammaWeightError, fgammaEmptyBinError = qcdPredictionHistos( data, plot, commonCut, True )

	egammaHist = getHists( data, commonCut, treeName="photonElectronTree" )
	egammaHist.Scale( fakeRate )

	fsrZ = getHists( [zgl,zgn], cut=commonCut )
	fsrW = getHists( [wg1,wg2], cut=commonCut )
	fsrTT = getHists( [tg], cut=commonCut )

	writeFinalTableLong( versionData, fakeRate, fakeRateStatError, fakeRateSysError, fakeRateSysErrorOwn,
		isrUncertaintyZ, isrUncertaintyW, isrUncertaintyT,
		dataHist, fgammaHist, fgammaWeightError, fgammaEmptyBinError,
		egammaHist, fsrZ, fsrW, fsrTT )

	# Add statistical uncertainties (ewk+yutaro)
	for bin in range( egammaHist.GetNbinsX()+2 ):
		egammaHist.SetBinError( bin, sqrt( egammaHist.GetBinError(bin)**2 +
			( egammaHist.GetBinContent(bin)/fakeRate*fakeRateStatError )**2 ) )

	# get all STATISTICAL uncertainties:
	statUncertHistStack = ROOT.THStack()
	statUncertHistStack.Add( fgammaHist )
	statUncertHistStack.Add( egammaHist )
	statUncertHistStack.Add( applyAddUncertainty( fsrZ.Clone(), 0 ) )
	statUncertHistStack.Add( applyAddUncertainty( fsrW.Clone(), 0 ) )
	statUncertHistStack.Add( applyAddUncertainty( fsrTT.Clone(), 0 ) )

	# and compure relative stat uncert
	statUncertHist = statUncertHistStack.GetStack().Last().Clone( randomName() )
	for bin in range( statUncertHist.GetNbinsX()+2 ):
		if statUncertHist.GetBinContent(bin):
			statUncertHist.SetBinError( bin, statUncertHist.GetBinError(bin) / statUncertHist.GetBinContent(bin) )
			statUncertHist.SetBinContent( bin, 1 )
	statUncertHist.SetFillStyle(3945)
	statUncertHist.SetMarkerSize(0)
	statUncertHist.SetFillColor(ROOT.kRed)


	# Now add up all uncertainties
	# QCD
	for bin in range( fgammaHist.GetNbinsX()+2 ):
		fgammaHist.SetBinError( bin, sqrt( fgammaHist.GetBinError(bin)**2
				+ fgammaWeightError.GetBinContent(bin)**2
				+ fgammaEmptyBinError.GetBinContent(bin)**2 ) )

	# EWK
	for bin in range( egammaHist.GetNbinsX()+2 ):
		egammaHist.SetBinError( bin, sqrt( egammaHist.GetBinError(bin)**2 +
			( egammaHist.GetBinContent(bin)/fakeRate )**2 * (fakeRateSysError**2+fakeRateSysErrorOwn**2 ) ) )

	# ISR
	fsrTT = applyAddUncertainty( fsrTT, isrUncertaintyT )
	fsrW = applyAddUncertainty( fsrW, isrUncertaintyW )
	fsrZ = applyAddUncertainty( fsrZ, isrUncertaintyZ )


	# Extract the signal and prettify the plots
	signal1 = getHists( ["slimW_1200_1120_375_V02.44_tree.root"], cut=commonCut )
	#signal2 = getHists( ["slimW_1400_1420_375_V01.00_tree.root"], cut=commonCut )
	#signal2.SetLineColor(ROOT.kMagenta-7)

	fsrZ.SetLineColor( ROOT.kRed-7 )
	fsrW.SetLineColor( ROOT.kRed-9 )
	fsrTT.SetLineColor( ROOT.kRed )
	signal1.SetLineColor( ROOT.kBlue )
	egammaHist.SetLineColor( 3 )

	for h in [fsrTT, signal1, fgammaHist, dataHist]:
		h.GetXaxis().SetTitle("#met#text{ [GeV]}")
		h.SetTitleSize(1./31.4485, "xy")
		h.SetLabelSize(1./31.4485, "xy")

	mh = Multihisto()
	mh.orderByIntegral = False
	mh.setMinimum(0.02)
	mh.leg.SetX1(0.75)
	mh.leg.SetX2(0.9)
	mh.leg.SetY1(0.65)
	mh.addHisto( fsrTT, "#gamma t#bar{t}", True )
	mh.addHisto( fsrW, "#gamma W", True )
	mh.addHisto( fsrZ, "#gamma Z", True )
	mh.addHisto( egammaHist, "e#rightarrow#gamma", True )
	mh.addHisto( fgammaHist, "QCD", True )

	mh.addHisto( dataHist, "Data", draw="pe" )
	mh.addHisto( signal1, "Signal", draw="hist" )
	#mh.addHisto( signal2, "1400/1420 Wino", draw="hist" )


	# draw stuff
	luminosity = 19.8
	infoText = ROOT.TLatex(0,.96, "#text{CMS Private Work#hspace{5cm}   }#SI{%s}{fb^{-1}}#, #sqrt{s}=#SI{8}{TeV}#, #geq1#gamma,#geq2#text{jets}"%luminosity )
	infoText.SetNDC()
	infoText.SetTextSize(1./31.4485)

	ROOT.gStyle.SetPaperSize(14.6,50.)
	ROOT.gStyle.SetPadTopMargin(0.05)
	ROOT.gStyle.SetPadRightMargin(0.02)
	ROOT.gStyle.SetPadLeftMargin(0.09)



	can = ROOT.TCanvas()
	mh.Draw()
	errorBand = mh.stack.GetStack().Last().Clone( randomName() )
	errorBand.SetFillStyle(3554)
	errorBand.SetFillColor(1)
	errorBand.SetMarkerSize(0)
	errorBand.Draw("e2 same")

	for ding in can.GetListOfPrimitives():
		if isinstance( ding, ROOT.TH1 ) or isinstance( ding, ROOT.THStack):
			ax = ding.GetYaxis()
			ax.SetTitleSize(1./31.4485)
			ax.SetLabelSize(1./31.4485)
			ax.SetTitleOffset(1.3)
			ax.SetLabelOffset(0)

	#printDeviations( dataHist, mh.stack.GetStack().Last() )

	from myRatio import Ratio
	r = Ratio( "Data / Bkg", dataHist, mh.stack.GetStack().Last() )
	r.draw(0,2)
	statUncertHist.Draw("same e2")
	r.ratio.Draw("same e")
	infoText.Draw()

	SaveAs( can, "finalDistributionData" )
	savePath = "/home/knut/master/documents/thesis/plots/finalDistributionData.tex"
	can.SaveAs( savePath )
	correctTiksPlot( savePath )



if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("--plot", nargs="+", default = ["met"] )
	opts = arguments.parse_args()

	if opts.plot == ["all"]:
		opts.plot = [ "met", "ht", "photons[0].ptJet()","Length$(jets.pt)", "Length$(photons.pt)"]

	for plot in opts.plot:
		finalDistributionData( plot )
