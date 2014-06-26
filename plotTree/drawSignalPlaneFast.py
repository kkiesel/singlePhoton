#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *
import re

Styles.tdrStyle2D()
ROOT.gStyle.SetOptLogz(0)

def interpolateEmptyBins( histo ):
	fillings = {}
	for x in range( 1, histo.GetNbinsX()+1 ):
		for y in range( 1, histo.GetNbinsY()+1 ):
			if histo.GetBinContent(x,y) == 0:
				values = [ histo.GetBinContent( x, y+1 ),
						histo.GetBinContent( x, y-1 ),
						histo.GetBinContent( x+1, y ),
						histo.GetBinContent( x-1, y ) ]
				values = [ v for v in values if v > 0 ]
				#if 'B' in histo.fileName and x == 16 and y == 15:
				#	values = [ histo.GetBinContent( x, y-1 ),
				#		histo.GetBinContent( x+1, y )]
				if len(values) == 4:
					fillings[(x,y)] = sum(values)/len(values)

	for (x,y), value in fillings.iteritems():
		histo.SetBinContent( x, y, value )

	return histo

def isSubset( l1, l2 ):
	#Returns true if all elements of l1 are in l2
	return len([ 1 for x in l1 if x in l2 ])==len(l1)

class Scan:
	def __init__( self ):
		self.points = {}
		self.functions = []

	def addHisto( self, x, y, name, histo ):
		if (x,y) not in self.points:
			self.points[(x,y)] = {}
		self.points[(x,y)][ name ] = histo

	def getGrid( self ):
		keys = self.points.keys()
		xList = sorted(set(zip(*keys)[0]))
		yList = sorted(set(zip(*keys)[1]))
		stepSizeX = xList[1] - xList[0]
		stepSizeY = yList[1] - yList[0]

		if xList != range( xList[0], xList[-1]+1, stepSizeX ):
			print "Error: no equidistant x-binning"
			print xList
			print "Assume T5wg sample and set bins manually"
			xList = range( xList[0], xList[-1]+1, 50 )
		if yList != range( yList[0], yList[-1]+1, stepSizeY ):
			print "Error: no equidistant y-binning"
			print yList

		defaultTitle = ""
		if len(xList) == 17 and len(yList) == 17:
			defaultTitle = ";m_{#tilde{q}} [GeV];m_{#tilde{g}} [GeV]"
		if len(xList) == 15 and len(yList) in [30,61]:
			defaultTitle = ";m_{#tilde{g}} [GeV];m_{#tilde{#chi}^{0}_{2}} [GeV]"
		self.defaultHisto = ROOT.TH2F( randomName(), defaultTitle, len(xList), xList[0]-0.5*stepSizeX, xList[-1]+0.5*stepSizeX, \
			len(yList), yList[0]-0.5*stepSizeY, yList[-1]+0.5*stepSizeY )

	def addFunction( self, function, name, histoTitle ):
		histo = self.defaultHisto.Clone( name )
		histo.GetZaxis().SetTitle( histoTitle )
		self.functions.append( (function, histo) )

	def fillFunctions( self ):
		for coordinate, histoDict in self.points.iteritems():
			bin = self.defaultHisto.FindBin( *coordinate )
			for function, histo in self.functions:
				histo.SetBinContent( bin, function( histoDict ) )


def getSignalHistosFromFile( filename ):
	# contains 'string', 'int', 'int', 'ROOT.TH1F'
	f = ROOT.TFile( filename )
	ROOT.gROOT.cd()
	histTuples = []
	for item in f.GetListOfKeys():
		name = item.GetName()
		if name in [ "photonTree", "photonElectronTree", "photonJetTree", "metFilters" ]: continue
		if not re.match( ".*\d+.*", name ): continue # if there is no digit in the name
		m = re.match( "([a-zA-Z]+)(\d+)_(\d+)", name )
		if not m:
			print "Could not match histogram name", name
		else: # put histogram in dictionary
			name, x, y = m.groups()
			x = int(x)
			y = int(y)
			histTuples.append( (name, x, y, item.ReadObj().Clone( name ) ) )
	return histTuples


############## functions used for creating histograms ###################################

def nGen( histoDict ):
	if "nGen" not in histoDict.keys(): return 0
	return histoDict["nGen"].GetBinContent(1)/1e4

def acceptance( histoDict, minMet=100, maxMet=-1 ):
	if "gMet" not in histoDict.keys(): return 0

	minBin = histoDict["gMet"].FindBin( minMet ) # corresponds to met cut
	if maxMet != -1: maxMet = histoDict["gMet"].FindBin( minMet )

	sel = histoDict["gMet"].Integral(minBin, maxMet )
	nGen = histoDict["nGen"].GetBinContent(1)

	return 100. * sel / nGen if nGen else 0

def signalContamination( histoDict ):
	if not isSubset( ["fMet", "eMet"], histoDict.keys() ): return 0
	minBin = histoDict["gMet"].FindBin( 100 ) # corresponds to met cut
	return 100.*addHistos([ histoDict["fMet"], histoDict["eMet"] ]).Integral(minBin, -1) / histoDict["gMet"].Integral(minBin, -1) if histoDict["gMet"].Integral(minBin, -1) else 0

def meanHt( histoDict ):
	if "gHt" not in histoDict.keys(): return 0
	return histoDict["gHt"].GetMean()

def meanPt( histoDict ):
	if "gPt" not in histoDict.keys(): return 0
	return histoDict["gPt"].GetMean()

def meanNjets( histoDict ):
	if "gNJets" not in histoDict.keys(): return 0
	return histoDict["gNJets"].GetMean()

def meanMet( histoDict ):
	if "gMet" not in histoDict.keys(): return 0
	return histoDict["gMet"].GetMean()

def jetScale( histoDict ):
	if not isSubset( ["gMetJecDown", "gMetJecUp"], histoDict.keys()): return 0

	#up = histoDict["gMetJecUp"].Integral()
	down = histoDict["gMetJecDown"].Integral()
	normal = histoDict["gMet"].Integral()
	#return 100.*(up-down)/(2*normal) if normal else 0
	return 100.*(normal-down)/normal if normal else 0

def puUncert( histoDict ):
	if not isSubset( ["gMetPuDown", "gMetPuUp" ], histoDict.keys()): return 0
	minMetValue = 100
	if isinstance( histoDict["gMetPuDown"], ROOT.TH1 ):
		minMetValueBin = histoDict["gMetPuDown"].FindBin( minMetValue )
		up = histoDict["gMetPuUp"].Integral( minMetValueBin, -1 )
		down = histoDict["gMetPuDown"].Integral( minMetValueBin, -1 )
		normal = histoDict["gMet"].Integral( minMetValueBin, -1 )
		return 100.*abs(up-down)/(2*normal) if normal else 0
	else:
		return 0

def statUncert( histoDict ):
	if "gMet" not in histoDict.keys(): return 0

	minMetValue = 350
	if isinstance( histoDict["gMet"], ROOT.TH1 ):
		bin = histoDict["gMet"].FindBin(350)
		i, e = integralAndError( histoDict["gMet"], bin, -1 )
		return 100*e/i if i else 0
	return 0

############## functions used for creating histograms end ################################





xSecPath = "Spectra_gsq_%s_8TeV.xsec"
pdfPath = "Spectra_gsq_%s_phad_pdfuncert.dat"


def xsection( self, file_ ):
	xsec = readSignalXSection( file_ )
	for x in self.xList:
		for y in self.yList:
			self.histo.SetBinContent( self.histo.FindBin(x,y), xsec[(x,y)][0] )
	self.histo.GetZaxis().SetTitle( "#sigma #text{ [pb]}" )
	self.draw( "xsectionSignal" )

def pdfUncertainty( self, file_ ):
	pdf = readSignalPdfUncertainty( file_ )
	for x in self.xList:
		for y in self.yList:
			xsec, acc = pdf[(x,y)]
			self.histo.SetBinContent( self.histo.FindBin(x,y), sqrt(xsec**2+acc**2) )
	self.histo.GetZaxis().SetTitle( "pdf uncert. [%]" )
	self.draw( "pdfUncertainty" )

def xsectionUncertainty( self, file_ ):
	xsec = readSignalXSection( file_ )
	for x in self.xList:
		for y in self.yList:
			xSec = xsec[(x,y)][0]
			uncertUp = xsec[(x,y)][1]
			uncertDown = xsec[(x,y)][2]
			self.histo.SetBinContent( self.histo.FindBin(x,y), (uncertUp+uncertDown)/2/xSec*100 )
	self.histo.GetZaxis().SetTitle( "#text{rel. }#sigma#text{ uncert.}#text{ [%]}" )
	self.draw( "xsectionSignalUncertaintyRel" )

def readSModelXsection( filename ):
	import re
	f = open( filename )

	xSections = {}

	for line in f.readlines():
		if line.startswith("#"): continue
		m, xsec, uncert = line.split(" ")
		xSections[int(m)] = (float(xsec), float(uncert) )
	f.close()

	return xSections


def sModelXsection( histogram, filename ):
	histogram.SetZTitle( "cross section [pb]" )
	uncert = histogram.Clone( "xSectionUncert" )
	uncert.SetZTitle( "cross section uncertainty [%]" )
	xSections = readSModelXsection( filename )
	for xBin in range(histogram.GetNbinsX()+2):
		binxSec = xSections[ int( histogram.GetXaxis().GetBinCenter(xBin) ) ]
		for yBin in range(histogram.GetNbinsY()+2):
			if histogram.GetYaxis().GetBinCenter(yBin) > histogram.GetXaxis().GetBinCenter(xBin): continue
			histogram.SetBinContent( xBin, yBin, binxSec[0] )
			uncert.SetBinContent( xBin, yBin, binxSec[1] )

	return [histogram, uncert]


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	opts = arguments.parse_args()

	for filename in opts.filenames:

		scanname = ""
		if "T5gg_" in filename:
			scanname = "T5gg"
		if "T5wg_" in filename:
			scanname = "T5wg"
		if "W_V03." in filename:
			scanname = "Wino"
		if "B_V03." in filename:
			scanname = "Bino"
		if not scanname: print "Cannot determine scan name"

		histTuples = getSignalHistosFromFile( filename )

		scan = Scan()
		for name, x, y, hist in histTuples:
			scan.addHisto( x, y, name, hist )
		scan.getGrid()

		scan.addFunction( nGen, "nGen", "Generated events [10^{4}]" )
		scan.addFunction( acceptance, "acceptance", "Acceptance [%]" )
		scan.addFunction( jetScale, "jes", "jes uncert. [%]" )
		scan.addFunction( signalContamination, "signalContamination", "bkg prediction from signal [%]" )
		scan.addFunction( meanHt, "ht", "#LTH_{T}#GT [GeV]" )
		scan.addFunction( meanPt, "pt", "#LTp_{T}#GT [GeV]" )
		scan.addFunction( meanNjets, "nJet", "#LTn_{Jets}#GT" )
		scan.addFunction( meanMet, "met", "#LT#slash{E}_{T}#GT [GeV]")
		scan.addFunction( puUncert, "pu", "pileup uncertainty [%]" )
		scan.addFunction( statUncert, "statUncert", "stat. uncert. [%]" )

		scan.fillFunctions()
		histos = list(zip( *scan.functions )[1])

		if scanname in ["T5gg", "T5wg" ]: histos.extend( sModelXsection( scan.defaultHisto.Clone( "xSection" ), "../../infos/simplifiedModel.xsec" ) )

		for histo in histos:
			can = ROOT.TCanvas()
			can.cd()

			#histo = interpolateEmptyBins( histo )
			name = histo.GetName()
			if name == "ht":
				histo.SetMinimum( 500 )
			elif name == "nJet":
				histo.SetMinimum( 2 )
			elif name == "pt":
				histo.SetMinimum( 110 )
			elif name == "met":
				histo.SetMinimum( 100 )
				histo.SetMaximum( 400 )
			elif name == "xSection":
				can.SetLogz()


			gr2D = ROOT.TGraph2D( histo )
			gr2D.SetMinimum( histo.GetMinimum() )
			gr2D.SetMaximum( histo.GetMaximum() )
			try:
				gr2D.SetTitle( "%s;%s;%s;%s"%(histo.GetTitle(),histo.GetXaxis().GetTitle(),histo.GetYaxis().GetTitle(),histo.GetZaxis().GetTitle()) )
				histo.GetZaxis().SetTitleOffset( 0.85 )
				gr2D.GetZaxis().SetTitleOffset( 0.85 )
			except:
				pass
			gr2D.SetNpx(100)
			gr2D.SetNpy(100)
			gr2D.Draw("colz")
			#histo.Draw("colz")

			scanText = ""
			if "T5gg_" in filename:
				scanText = "- pp#rightarrow#tilde{g}#tilde{g}#rightarrow#tilde{#chi}#tilde{#chi}#gamma#gamma,"
			if "T5wg_" in filename:
				scanText = "- pp#rightarrow#tilde{g}#tilde{g}#rightarrow#tilde{#chi}#tilde{#chi}W#gamma,"


			text = "CMS Private Work %s 19.7fb^{-1} 8TeV #geq1#gamma#geq2jets"%scanText
			if name in ["xSection", "xSectionUncert"]:
				text = "CMS Private Work %s 19.7fb^{-1} 8TeV #geq1#gamma#geq2jets"%"pp#rightarrow#tilde{g}#tilde{g}"

			info = ROOT.TLatex(.01,.96, text )
			info.SetNDC()
			info.Draw()

			ROOT.gPad.SaveAs( "plots/%s_%s.pdf"%(scanname,histo.GetName()) )

#todo xsectiond, pdfuncert, xsectionuncert
#todo ngen correct

