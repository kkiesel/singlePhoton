#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *
import re

Styles.tdrStyle2D()
ROOT.gStyle.SetOptLogz(0)
ROOT.gStyle.SetPalette(1)

def nGenHisto( nGen ):
	h = ROOT.TH1F(randomName(), ";met;", 1, 0, 1 )
	h.SetBinContent( 1, nGen )
	return h

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
				if len(values) >= 3:
					fillings[(x,y)] = sum(values)/len(values)

	for (x,y), value in fillings.iteritems():
		histo.SetBinContent( x, y, value )

	return histo

def isSubset( l1, l2 ):
	#Returns true if all elements of l1 are in l2
	return len([ 1 for x in l1 if x in l2 ])==len(l1)

def histoToList( histo, uncert=False ):
	import array
	binning = [ 100, 120, 160, 200, 270, 350 ]
	metBinning = array.array( "d", binning )
	h = histo.Rebin( len(metBinning)-1, randomName(), metBinning )
	out = []
	for bin in range( 1, h.GetNbinsX()+2 ):
		if uncert:
			out.append( h.GetBinError(bin) )
		else:
			out.append( h.GetBinContent(bin) )
	return out

def histoDiffToList( histo1, histo2 ):
	import array
	binning = [ 100, 120, 160, 200, 270, 350 ]
	metBinning = array.array( "d", binning )
	h1 = histo1.Rebin( len(metBinning)-1, randomName(), metBinning )
	h2 = histo2.Rebin( len(metBinning)-1, randomName(), metBinning )
	out = []
	h1.Add( h2, -1 )
	h1.Scale(0.5) # (h1 - h2)/2
	for bin in range( 1, h1.GetNbinsX()+2 ):
		out.append( abs(h1.GetBinContent(bin)) )
	return out

def listToString( list ):
	return " ".join([str(i) for i in list ] )

class Scan:
	def __init__( self ):
		self.points = {}
		self.functions = []
		self.nGen = -1
		self.binning = [ 100, 120, 160, 200, 270, 350, -1 ]
		self.results = []

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

	def addFunction( self, function, name, histoTitle, splitBins=False ):
		histo = self.defaultHisto.Clone( name )
		histo.GetZaxis().SetTitle( histoTitle )
		self.functions.append( (function, histo, None) )
		if splitBins:
			for i, minMet in enumerate( self.binning[0:-1] ):
				maxMet = self.binning[i+1]
				histo = self.defaultHisto.Clone( name+str(i) )
				histo.GetZaxis().SetTitle( histoTitle )
				self.functions.append( (function, histo, (minMet,maxMet) ) )


	def fillFunctions( self ):
		for coordinate, histoDict in self.points.iteritems():
			#if "nGen" not in histoDict:
			if nGen != -1:
				histoDict["nGen"] = nGenHisto( self.nGen )

			bin = self.defaultHisto.FindBin( *coordinate )
			for function, histo, options in self.functions:
				if options:
					histo.SetBinContent( bin, function( histoDict, *options ) )
				else:
					histo.SetBinContent( bin, function( histoDict ) )

	def fillResults( self, filename ):

		infos = []
		counter = 0
		import time
		signalCardString = """
# This file contains event yield and uncertainties for signal points.
# Produced by:
#    Knut Kiesel
# on:
#    %s
# with input file:
#    %s
# There will be a common section for all signal points, folllowed by the information
# for each signal point. The numbering of the signal points is arbitrary.
# For bin edges, the ROOT conversion is used (lower bin edge included, upper bin
# edge excluded).
"""%( time.strftime("%Y-%m-%d %H:%M:%S"), filename )

		for coordinate, histoDict in self.points.iteritems():
			info = {}
			info["m1"] = coordinate[0]
			info["m2"] = coordinate[1]
			info["nGen"] = int(histoDict["nGen"].GetBinContent(1))
			info["signal"] = histoToList( histoDict["gMet"] )
			info["signalStat"] = histoToList( histoDict["gMet"], True )
			info["signalSystJes"] = histoDiffToList( histoDict["gMetJesUp"], histoDict["gMetJesDown"] ) if "gMetJesUp" in histoDict else []
			info["signalSystPu"] = histoDiffToList( histoDict["gMetPuUp"], histoDict["gMetPuDown"] ) if "gMetPuUp" in histoDict else []
			info["signalEWK"] = histoToList( histoDict["eMet"] )
			info["signalEWKStat"] = histoToList( histoDict["eMet"], True )
			info["signalEWKSyst"] = [0.11*i for i in histoToList( histoDict["eMet"] ) ] # 11% uncertainy
			info["signalQCD"] = histoToList( histoDict["fMet"] )
			info["signalQCDStat"] = histoToList( histoDict["fMet"], True )
			info["signalQCDSyst"] = histoToList( histoDict["fMetError"] )
			infos.append(info)

			signalsystUncert = [ sqrt(i**2+j**2) for i,j in zip( info["signalSystJes"], info["signalSystPu"] ) ]

			signalCardString += "Point %s gluino mass = %s\n"%(counter, info["m1"])
			signalCardString += "Point %s nlsp mass = %s\n"%(counter, info["m2"])
			signalCardString += "Point %s generated events = %s\n"%(counter, info["nGen"])
			signalCardString += "Point %s number of signal events in bins = %s 0\n"%(counter, listToString( info["signal"] ) )
			signalCardString += "Point %s statistical error of signal events in bins = %s 0\n"%(counter, listToString( info["signalStat"] ) )
			signalCardString += "Point %s systematical error of signal events in bins = %s 0\n"%(counter, listToString( signalsystUncert ) )
			signalCardString += "Point %s EWK prediction = %s 0\n"%(counter, listToString( info["signalEWK"] ) )
			signalCardString += "Point %s QCD prediction = %s 0\n"%(counter, listToString( info["signalQCD"] ) )
			signalCardString += "\n"

			counter += 1


		outputFileName = "eventYield%s-%s.txt"%(filename, time.strftime("%Y-%m-%d"))
		signalCardFile = open( outputFileName, "w")
		signalCardFile.write( signalCardString )
		signalCardFile.close()
		print "Wrote to %s"%outputFileName




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

def signalContamination( histoDict, minMet=100, maxMet=-1 ):
	if not isSubset( ["fMet", "eMet"], histoDict.keys() ): return 0

	minBin = histoDict["gMet"].FindBin( minMet )
	if maxMet != -1: maxMet = histoDict["gMet"].FindBin( minMet )

	minBin = histoDict["gMet"].FindBin( minBin, maxMet )
	return 100.*addHistos([ histoDict["fMet"], histoDict["eMet"] ]).Integral(minBin, maxMet) / histoDict["gMet"].Integral(minBin, maxMet) if histoDict["gMet"].Integral(minBin, maxMet) else 0

def meanHt( histoDict ):
	if "gHt" not in histoDict.keys(): return 0
	return histoDict["gHt"].GetMean()

def meanPt( histoDict ):
	if "gPt" not in histoDict.keys(): return 0
	return histoDict["gPt"].GetMean()

def meanNjets( histoDict ):
	if "gNjets" not in histoDict.keys(): return 0
	return histoDict["gNjets"].GetMean()

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

############## functions used for creating histograms end #####################
############## begin of declaration of other functions ########################

def gmsbPdfUncert( histogram, filename ):
	histogram = histogram.Clone( "pdfUncert" )
	histogram.SetZTitle( "pdf uncert. [%]" )
	pdfUncert = readSignalPdfUncertainty( filename )
	for multiBin, info in pdfUncert.iteritems():
		bin = histogram.FindBin( multiBin[0], multiBin[1] )
		histogram.SetBinContent( bin, sqrt(info[0]**2+info[1]**2) )
	return [histogram]

def gmsbXsection( histogram, filename ):
	histogram = histogram.Clone( "xSection" )
	histogram.SetZTitle( "cross section [pb]" )
	uncert = histogram.Clone( "xSectionUncert" )
	uncert.SetZTitle( "cross section uncertainty [%]" )
	xSections = readSignalXSection( filename )
	for multiBin, info in xSections.iteritems():
		bin = histogram.FindBin( multiBin[0], multiBin[1] )
		histogram.SetBinContent( bin, info[0] )
		uncert.SetBinContent( bin, 100*abs(info[1]-info[2])/info[0] )

	return [histogram, uncert]

def sModelXsection( histogram, filename ):
	histogram = histogram.Clone( "xSection" )
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

############## end of declaration of other functions ##########################

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

		if "W_gsq_" in filename:
			scanname = "Wino"
		if "B_gsq_" in filename:
			scanname = "Bino"

		if not scanname: print "Cannot determine scan name"

		histTuples = getSignalHistosFromFile( filename )

		scan = Scan()
		for name, x, y, hist in histTuples:
			scan.addHisto( x, y, name, hist )
		scan.getGrid()

		scan.fillResults( scanname )

		if scanname == "Bino": scan.nGen = 10000
		if scanname == "Wino": scan.nGen = 60000

		scan.addFunction( nGen, "nGen", "Generated events [10^{4}]" )
		scan.addFunction( acceptance, "acceptance", "Acceptance [%]", True )
		scan.addFunction( jetScale, "jes", "jes uncert. [%]" )
		scan.addFunction( signalContamination, "signalContamination", "rel. bkg. pred. from signal [%]", True )
		scan.addFunction( meanHt, "ht", "#LTH_{T}#GT [GeV]" )
		scan.addFunction( meanPt, "pt", "#LTp_{T}#GT [GeV]" )
		scan.addFunction( meanNjets, "nJet", "#LTn_{Jets}#GT" )
		scan.addFunction( meanMet, "met", "#LT#slash{E}_{T}#GT [GeV]")
		scan.addFunction( puUncert, "pu", "pileup uncertainty [%]" )
		scan.addFunction( statUncert, "statUncert", "stat. uncert. [%]" )

		scan.fillFunctions()
		histos = list(zip( *scan.functions )[1])

		if "T5" in scanname:
			histos.extend( sModelXsection( scan.defaultHisto, "../../infos/simplifiedModel.xsec" ) )

		if scanname[1:] == "ino" : # Bino or Wino
			scanAbbr = scanname[0] # B or W
			histos.extend( gmsbXsection( scan.defaultHisto, "../../infos/Spectra_gsq_%s_8TeV.xsec"%scanAbbr ) )
			histos.extend( gmsbPdfUncert( scan.defaultHisto, "../../infos/Spectra_gsq_%s_phad_pdfuncert.dat"%scanAbbr) )

		for histo in histos:
			can = ROOT.TCanvas()
			can.cd()

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
			if histo.GetEffectiveEntries():
				gr2D.SetMinimum( histo.GetMinimum() )
				gr2D.SetMaximum( histo.GetMaximum() )

			try:
				gr2D.SetTitle( "%s;%s;%s;%s"%(histo.GetTitle(),histo.GetXaxis().GetTitle(),histo.GetYaxis().GetTitle(),histo.GetZaxis().GetTitle()) )
				histo.GetZaxis().SetTitleOffset( 0.85 )
				if histo.GetEffectiveEntries():
					gr2D.GetZaxis().SetTitleOffset( 0.85 )
			except:
				pass
			gr2D.SetNpx(100)
			gr2D.SetNpy(100)
			gr2D.Draw("colz")
			if True:
				histo = interpolateEmptyBins( histo )
				histo.Draw("colz")

			scanText = ""
			if "T5gg_" in filename:
				scanText = "pp#rightarrow#tilde{g}#tilde{g}#rightarrow#tilde{#chi}#tilde{#chi}#gamma#gamma"
			if "T5wg_" in filename:
				scanText = "pp#rightarrow#tilde{g}#tilde{g}#rightarrow#tilde{#chi}#tilde{#chi}W#gamma"

			if scanname == "Wino":
				scanText = "Wino #tilde{#chi}_{1}"
			if scanname == "Bino":
				scanText = "Bino #tilde{#chi}_{1}^{0}"



			minMet = 100
			maxMet = -1

			if name[-1].isdigit():
				bin = int(name[-1])
				minMet = scan.binning[bin]
				maxMet = scan.binning[bin+1]

			cutText = "#geq1#gamma#geq2jets"
			if name in  [ "xSection", "xSectionUncert", "pdfUncert" ]:
				cutText = ""
				minMet = -1
				maxMet = -1

			if maxMet < 0 and minMet < 0:
				pass
			elif maxMet < 0:
				cutText +=" #slash{E}_{T}#geq%sGeV"%minMet
			else:
				cutText += " %s#leq#slash{E}_{T}(GeV)<%s"%(minMet,maxMet)


			info = ROOT.TLatex(.0,.96, "CMS Simulation - "+scanText+" - 8TeV" )
			info.SetNDC()
			info.SetTextSize( histo.GetLabelSize() )
			info.SetTextFont( histo.GetLabelFont() )
			info.Draw()
			info2 = info.Clone()
			info2.SetTextAlign(31) # right align
			info2.SetTextSize( info2.GetTextSize() * .9 )
			info2.DrawLatex( 0.98, 0.96, cutText )
			info2.Draw()

			ROOT.gPad.SaveAs( "plots/%s_%s.pdf"%(scanname,histo.GetName()) )
			#ROOT.gPad.SaveAs( "plots/%s_%s.png"%(scanname,histo.GetName()) )
			#ROOT.gPad.SaveAs( "plots/%s_%s.C"%(scanname,histo.GetName()) )

#todo ngen correct

