#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *
st = Styles.tdrStyle2D()
st.SetOptLogz(0)

class SignalPlane:
	def __init__( self, fileName ):
		self.histo = ROOT.TH2F( randomName(), ";m_{#tilde{q}}  [GeV];m_{#tilde{g}}  [GeV]", 17, 350, 2050, 17, 370, 2070 )
		self.histo.GetZaxis().SetTitleOffset( 0.85 )
		self.xList = range( 400, 2001, 100 )
		self.yList = range( 420, 2021, 100 )
		self.fileName = fileName

	def xsection( self, file_ ):
		xsec = readSignalXSection( file_ )
		for x in self.xList:
			for y in self.yList:
				self.histo.SetBinContent( self.histo.FindBin(x,y), xsec[(x,y)][0] )
		self.histo.GetZaxis().SetTitle( "#sigma  [pb]" )
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
		self.histo.GetZaxis().SetTitle( "rel. #sigma uncert. [%]" )
		self.draw( "xsectionSignalUncertaintyRel" )

	def xsecLimit( self, abbr ):
		histo = wHistoLimit if abbr == "W" else bHistoLimit
		for x in self.xList:
			for y in self.yList:
				bin = histo.FindBin(x,y)
				self.histo.SetBinContent( self.histo.FindBin(x,y), histo.GetBinContent(bin) )
		self.histo.GetZaxis().SetTitle( "expected #sigma limit [pb]" )
		self.draw( "excludedXsecLimit" )


	def fill( self, function, zlabel, saveName ):
		for x in self.xList:
			for y in self.yList:

				if "W_gsq_" in self.fileName:
					if x == 500 and y == 520 or x == 1000 and y == 1520:
						self.histo.SetBinContent( self.histo.FindBin(x,y), 0 )
						continue # these points are not available for wino sample

				if "B_gsq_" in self.fileName:
					if x == 1900 and y == 1820:
						self.histo.SetBinContent( self.histo.FindBin(x,y), 0 )
						continue

				histoNames = ["gMet", "eMet", "fMet", "fMetError", "gHt", "gPt", "gNjets", "gMetPuUp", "gMetPuDown", "gMetJecUp", "gMetJecDown" ]

				histos = {}
				for histoName in histoNames:
					histos[ histoName ] = readHisto( self.fileName, "%s%s_%s"%(
							histoName, x, y ) )
				#histos["eMet"] = applyFakeRateEWK( histos["eMet"] )

				# Save the number of generated events. In fact this is not a histo,
				# but it's easier to do so.
				if "acceptance" in saveName:
					histos["nGen"] = 60000 if "W_gsq" in self.fileName else 10000

				self.histo.SetBinContent( self.histo.FindBin(x,y), function( histos ) )
		self.histo.GetZaxis().SetTitle( zlabel )
		self.draw( saveName )

	def draw( self, saveName ):

		fillEmptyBins = True
		if fillEmptyBins:
			fillings = {}
			for x in range( 1, self.histo.GetNbinsX()+1 ):
				for y in range( 1, self.histo.GetNbinsY()+1 ):
					if self.histo.GetBinContent(x,y) == 0:
						values = [ self.histo.GetBinContent( x, y+1 ),
								self.histo.GetBinContent( x, y-1 ),
								self.histo.GetBinContent( x+1, y ),
								self.histo.GetBinContent( x-1, y ) ]
						values = [ v for v in values if v > 0 ]
						if 'B' in self.fileName and x == 16 and y == 15:
							values = [ self.histo.GetBinContent( x, y-1 ),
								self.histo.GetBinContent( x+1, y )]
						if len(values):
							fillings[(x,y)] = sum(values)/len(values)

			for (x,y), value in fillings.iteritems():
				self.histo.SetBinContent( x, y, value )
				pass

		ROOT.gStyle.SetPadRightMargin(0.2)
		can = ROOT.TCanvas()
		can.cd()
		if saveName in [ "xsectionSignal", "excludedXsecLimit" ]:
			can.SetLogz(1)
		else:
			can.SetLogz(0)
		self.histo.GetXaxis().SetNdivisions(5,5,0,True)
		# When the first Palette object is created (Draw("colz")), it will inherit
		# it's title from the histogram. If the histogram is changing, the palette
		# won't recognize.
		palette = self.histo.GetListOfFunctions().FindObject("palette")
		if palette:
			palette.GetAxis().SetTitle( self.histo.GetZaxis().GetTitle() )
		self.histo.GetZaxis().SetTitleOffset(1.2)

		# set range
		minimum = 1e6
		maximum = 0
		if saveName == "acceptance":
			minimum = 0
		if saveName == "signalPlaneHt":
			minimum = 500
			maximum = 1600*1.05
		if saveName == "signalPlanePt":
			minimum = 110
			maximum = 600*1.05
		if saveName == "signalPlaneMet":
			minimum = 100
			maximum = 350*1.05
		if saveName == "signalPlaneNjets":
			minimum = 2
			maximum = 6*1.05
		if saveName == "signalContaminationPlane":
			minimum = 0
			maximum = 3.6*1.05
		if minimum < 1e5:
			self.histo.SetMinimum( minimum )
		if maximum > 0:
			self.histo.SetMaximum( maximum )

		self.histo.Draw("colz")

		#"#text{CMS Private Work  -  }#SI{19.8}{fb^{-1}}#, #sqrt{s}=#SI{8}{TeV}#, #geq1#ggamma,#geq2#text{jets}#,, #met<#SI{100}{GeV}" )
		text = "CMS Private Work - "
		particleString = "Bino" if self.fileName[0] == "B" else "Wino"

		text = "%s-like #tilde{#chi}^{0}_{1}"%particleString

		if "xsectionSignal" in saveName:
			text += " 8TeV"
		elif "excludedXsecLimit" in saveName:
			text += "19.7fb^{-1} #geq1#gamma,#geq2jets"
		else:
			text += " 8TeV #geq1#gamma,#geq2jets"

		info = PlotCaption(treeName="")
		info = ROOT.TLatex(.02,.97, text )
		info.SetTextSize( info.GetTextSize()*0.9 )
		info.SetNDC()
		info.Draw()

		SaveAs( can, saveName+'_'+self.fileName[0] )


def signalContamination( h ):
	minBin = h["gMet"].FindBin( 100 ) # corresponds to met cut
	return 100.*addHistos([ h["fMet"], h["eMet"] ]).Integral(minBin, -1) / h["gMet"].Integral(minBin, -1) if h["gMet"].Integral(minBin, -1) else 0

def acceptance( h ):
	minBin = h["gMet"].FindBin( 100 ) # corresponds to met cut
	return 1.*h["gMet"].Integral(minBin, -1 )/h["nGen"]*100

def acceptanceLast( h ):
	minBin = h["gMet"].FindBin( 350 ) # corresponds to last met bin
	return 1.*h["gMet"].Integral(minBin, -1 )/h["nGen"]*100

def meanHt( h ):
	return h["gHt"].GetMean()

def meanPt( h ):
	return h["gPt"].GetMean()

def meanNjets( h ):
	return h["gNjets"].GetMean()

def meanMet( h ):
	return h["gMet"].GetMean()

def jetScale( h ):
	minMetValue = 100
	if isinstance( h["gMetJecDown"], ROOT.TH1 ):
		minMetValueBin = h["gMetJecDown"].FindBin( minMetValue )
		up = h["gMetJecUp"].Integral( minMetValueBin, -1 )
		down = h["gMetJecDown"].Integral( minMetValueBin, -1 )
		normal = h["gMet"].Integral( minMetValueBin, -1 )
		return 100.*(up-down)/(2*normal) if normal else 0
	else:
		return 0

def puUncert( h ):
	minMetValue = 100
	if isinstance( h["gMetPuDown"], ROOT.TH1 ):
		minMetValueBin = h["gMetPuDown"].FindBin( minMetValue )
		up = h["gMetPuUp"].Integral( minMetValueBin, -1 )
		down = h["gMetPuDown"].Integral( minMetValueBin, -1 )
		normal = h["gMet"].Integral( minMetValueBin, -1 )
		return 100.*abs(up-down)/(2*normal) if normal else 0
	else:
		return 0

def statUncert( h ):
	minMetValue = 350
	if isinstance( h["gMet"], ROOT.TH1 ):
		bin = h["gMet"].FindBin(350)
		i, e = integralAndError( h["gMet"], bin, -1 )
		return 100*e/i if i else 0
	return 0

sp = SignalPlane( "W_V03.02.root" )
sp.fill( jetScale, "jet energy scale uncertainty [%]", "jesUncert" )
sp = SignalPlane( "W_V03.02.root" )
sp.fill( puUncert, "pile up uncertainty [%]", "puUncert" )
sp = SignalPlane( "W_V03.02.root" )
sp.fill( statUncert, "stat. uncert. [%] #slash{E}_{T}#geq350GeV", "statUncert" )

import sys
sys.exit()

#filePath = "%s_gsq_V02.45.root"
filePath = "%s_V03.01.root"
xSecPath = "Spectra_gsq_%s_8TeV.xsec"
pdfPath = "Spectra_gsq_%s_phad_pdfuncert.dat"

#for scanAbbr in [ "B", "W" ]:
#	sp = SignalPlane( filePath%scanAbbr )
#	sp.pdfUncertainty( pdfPath%scanAbbr )
#	sp.xsection( xSecPath%scanAbbr )
#	sp.xsectionUncertainty( xSecPath%scanAbbr )
#	sp.xsecLimit( scanAbbr )

#for scanAbbr in [ "B", "W" ]:
#	sp = SignalPlane( filePath%scanAbbr )
#	sp.fill( acceptanceLast, "Acceptance last Bin [%]", "acceptanceLastBin" )

for scanAbbr in [ "B", "W" ]:
	sp = SignalPlane( filePath%scanAbbr )
	sp.fill( acceptance, "Acceptance [%]", "acceptance" )

for scanAbbr in [ "B", "W" ]:
	sp = SignalPlane( filePath%scanAbbr )
	sp.fill( signalContamination, "Signal Contamination [%]", "signalContaminationPlane" )

for scanAbbr in [ "B", "W" ]:
	sp = SignalPlane( filePath%scanAbbr )
	sp.fill( meanHt, "H_{T} [GeV]", "signalPlaneHt" )

for scanAbbr in [ "B", "W" ]:
	sp = SignalPlane( filePath%scanAbbr )
	sp.fill( meanPt, "p_{T^{*}} [GeV]", "signalPlanePt" )

for scanAbbr in [ "B", "W" ]:
	sp = SignalPlane( filePath%scanAbbr )
	sp.fill( meanNjets, "n_{jets}", "signalPlaneNjets" )

for scanAbbr in [ "B", "W" ]:
	sp = SignalPlane( filePath%scanAbbr )
	sp.fill( meanMet, "#met [GeV]", "signalPlaneMet" )

	#sp.fill( signalContaminationError, "#sigma_{signal contermination}", "signalConmatinationError" )
	#sp.fill( signalContaminationErrorRel, "#sigma_{signal cont.} / signal cont.", "signalConmatinationErrorRel" )
	#sp.fill( signalContaminationE, "signal contamination e#rightarrow#gamma", "signalContaminationPlaneE" )

