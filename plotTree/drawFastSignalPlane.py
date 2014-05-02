#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *

def signalContamination( h ):
	minBin = h["gMet"].FindBin( 100 ) # corresponds to met cut
	scaledEMet = applyFakeRateEWK( h["eMet"].Clone( randomName() ) )
	return 100.*addHistos([ h["fMet"], scaledEMet ]).Integral(minBin, -1) / h["gMet"].Integral(minBin, -1) if h["gMet"].Integral(minBin, -1) else 0

def acceptance( h ):
	minBin = h["gMet"].FindBin( 100 ) # corresponds to met cut
	return 1.*h["gMet"].Integral(minBin, -1 )/h["nGen"]*100

def meanHt( h ):
	return h["gHt"].GetMean()

def meanPt( h ):
	return h["gPt"].GetMean()

def meanNjets( h ):
	return h["gNjets"].GetMean()

def meanMet( h ):
	return h["gMet"].GetMean()

class SignalPlaneLooper:
	def __init__( self, fileName ):
		self.defaultHisto = ROOT.TH2F( randomName(), ";m_{#tilde{q}} #text{ [GeV]};m_{#tilde{g}} #text{ [GeV]}", 17, 350, 2050, 17, 370, 2070 )
		self.xList = range( 400, 2001, 100 )
		self.yList = range( 420, 2021, 100 )
		self.fileName = fileName
		self.funcDict = {}
		self.histDict = {}

	def iniHisto( self, abbr, function ):
		self.histDict[abbr] = self.defaulHisto.Clone( randomName() )
		self.funcDict[abbr] = function

	def loop( self ):

		histoNames = ["gMet", "eMet", "fMet", "fMetUp", "fMetDown", "gHt", "gPt", "gNjets" ]

		for x in self.xList:
			for y in self.yList:

				# read in histograms
				histos = {}
				for histoName in histoNames:
					try:
						histos[ histoName ] = readHisto( self.fileName, "%s%s_%s"%( histoName, x, y ) )
					except:
						print "Did not found histogram %s in file %s"%histoName, self.fileName
						continue

				# use functions
				for abbr, hist in self.histDict.iteritems():
					hist.SetBinContent( hist.FindBin(x,y), self.funcDict[abbr](histos) )

def fillEmptyBins( histo ):
	# fills empty bins with mean of the four neighbours
	fillings = {}
	for x in range( 1, histo.GetNbinsX()+1 ):
		for y in range( 1, histo.GetNbinsY()+1 ):
			if histo.GetBinContent(x,y) == 0:
				values = [ histo.GetBinContent( x, y+1 ),
						histo.GetBinContent( x, y-1 ),
						histo.GetBinContent( x+1, y ),
						histo.GetBinContent( x-1, y ) ]
				values = [ v for v in values if v > 0 ]

				if x == 16 and y == 15: # special case for bino, near diagonal
					values = [ histo.GetBinContent( x, y-1 ),
						histo.GetBinContent( x+1, y )]
				if len(values):
					fillings[(x,y)] = sum(values)/len(values)

	for (x,y), value in fillings.iteritems():
		histo.SetBinContent( x, y, value )

	return histo




filePath = "%s_gsq_V02.45.root"
xSecPath = "/home/knut/master/infos/Spectra_gsq_%s_8TeV.xsec"
pdfPath = "/home/knut/master/infos/Spectra_gsq_%s_phad_pdfuncert.dat"

for scanAbbr in [ "B", "W" ]:
	sp = SignalPlaneLooper( filePath%scanAbbr )
	sp.iniHisto( "acc", acceptance )
	sp.iniHisto( "ht", meanHt )
	sp.iniHisto( "met", meanMet )

	sp.loop()


st = Styles.tdrStyle2D()
st.SetOptLogz(0)
paperWidth = 14.65 #cm
paperWidth2 = 5.7
st.SetPaperSize(paperWidth2,50.)

