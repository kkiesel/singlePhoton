#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *
st = Styles.tdrStyle2D()
st.SetOptLogz(0)
paperWidth = 14.65 #cm
paperWidth2 = 5.7
st.SetPaperSize(paperWidth2,50.)


class SignalPlane:
	def __init__( self, fileName ):
		self.histo = ROOT.TH2F( randomName(), ";m_{#tilde{q}} #text{ [GeV]};m_{#tilde{g}} #text{ [GeV]}", 17, 350, 2050, 17, 370, 2070 )
		self.histo.GetZaxis().SetTitleOffset( 0.85 )
		self.xList = range( 400, 2001, 100 )
		self.yList = range( 420, 2021, 100 )
		self.fileName = fileName

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
		self.histo.GetZaxis().SetTitle( "pdf uncert." )
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

				histoNames = ["gMet", "eMet", "fMet", "fMetUp", "fMetDown", "gHt", "gPt", "gNjets" ]

				histos = {}
				for histoName in histoNames:
					histos[ histoName ] = readHisto( self.fileName, "%s%s_%s"%(
							histoName, x, y ) )
				histos["eMet"] = applyFakeRateEWK( histos["eMet"] )

				# Save the number of generated events. In fact this is not a histo,
				# but it's easier to do so.
				if "acceptance" in saveName:
					histos["nGen"] = 60000 if "W_gsq" in self.fileName else 10000

				self.histo.SetBinContent( self.histo.FindBin(x,y), function( histos ) )
		self.histo.GetZaxis().SetTitle( zlabel )
		self.draw( saveName )

	def draw( self, saveName ):
		ROOT.gStyle.SetPadRightMargin(0.2)
		can = ROOT.TCanvas()
		can.cd()
		if "xsectionSignal" == saveName:
			can.SetLogz(1)
		else:
			can.SetLogz(0)
		labelScale = 1.*0.49062
		titleScale = 1.*0.588744
		self.histo.SetTitleSize(0.06/titleScale, "XYZ")
		self.histo.SetLabelSize(0.05/labelScale, "XYZ")
		self.histo.SetLabelOffset(0.0001, "XYZ" )
		self.histo.SetTitleOffset(0.5,"z" )
		self.histo.GetXaxis().SetNdivisions(5,5,0,True)
		# When the first Palette object is created (Draw("colz")), it will inherit
		# it's title from the histogram. If the histogram is changing, the palette
		# won't recognize.
		palette = self.histo.GetListOfFunctions().FindObject("palette")
		if palette:
			palette.GetAxis().SetTitle( self.histo.GetZaxis().GetTitle() )
		self.histo.GetZaxis().SetTitleOffset(1.2)
		self.histo.Draw("colz")

		#"#text{CMS Private Work  -  }#SI{19.8}{fb^{-1}}#, #sqrt{s}=#SI{8}{TeV}#, #geq1#gamma,#geq2#text{jets}#,, #met<#SI{100}{GeV}" )
		text = "#text{CMS Private Work - }"
		particleString = "Bino" if self.fileName[0] == "B" else "Wino"

		text = "#text{%s-like }#tilde{#chi}^0_1"%particleString

		if "xsectionSignal" in saveName:
			text += "#, #sqrt{s}=#SI{8}{TeV}"
		else:
			text += "#, #sqrt{s}=#SI{8}{TeV}#, #geq1#gamma,#geq2#text{jets}"

		info = PlotCaption(treeName="")
		info = ROOT.TLatex(-.12,.98, text )
		info.SetNDC()
		info.SetTextSize(1/9.81241)
		info.Draw()

		SaveAs( can, saveName+'_'+self.fileName[0] )

		texFileName = "/home/knut/master/documents/thesis/plots/%s_%s.tex"%(saveName,self.fileName[0])
		can.SaveAs( texFileName )
		correctTiksPlot( texFileName )



def signalContamination( h ):
	minBin = h["gMet"].FindBin( 100 ) # corresponds to met cut
	return 100.*addHistos([ h["fMet"], h["eMet"] ]).Integral(minBin, -1) / h["gMet"].Integral(minBin, -1) if h["gMet"].Integral(minBin, -1) else 0

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

filePath = "%s_gsq_V02.45.root"
xSecPath = "/home/knut/master/infos/Spectra_gsq_%s_8TeV.xsec"
pdfPath = "/home/knut/master/infos/Spectra_gsq_%s_phad_pdfuncert.dat"

for scanAbbr in [ "B", "W" ]:
	sp = SignalPlane( filePath%scanAbbr )
	sp.pdfUncertainty( pdfPath%scanAbbr )
	sp.xsection( xSecPath%scanAbbr )
	sp.xsectionUncertainty( xSecPath%scanAbbr )
	sp.fill( acceptance, "Acceptance [%]", "acceptance" )
	sp.fill( signalContamination, "Signal Contamination [%]", "signalContaminationPlane" )
	sp.fill( meanHt, "H_{T}#text{ [GeV]}", "signalPlaneHt" )
	sp.fill( meanPt, "p_{T^{*}}#text{ [GeV]}", "signalPlanePt" )
	sp.fill( meanNjets, "n_{#text{jets}}", "signalPlaneNjets" )
	sp.fill( meanMet, "#met#text{ [GeV]}", "signalPlaneMet" )


	#sp.fill( signalContaminationError, "#sigma_{signal contermination}", "signalConmatinationError" )
	#sp.fill( signalContaminationErrorRel, "#sigma_{signal cont.} / signal cont.", "signalConmatinationErrorRel" )
	#sp.fill( signalContaminationE, "signal contamination e#rightarrow#gamma", "signalContaminationPlaneE" )

