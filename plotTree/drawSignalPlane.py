#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *
Styles.tdrStyle2D()

class SignalPlane:
	def __init__( self, fileName ):
		self.histo = ROOT.TH2F( randomName(), ";M_{#tilde{q}} [GeV];M_{#tilde{g}} [GeV]", 17, 350, 2050, 17, 370, 2070 )
		self.histo.GetZaxis().SetTitleOffset( 0.85 )
		self.xList = range( 400, 2001, 100 )
		self.yList = range( 420, 2021, 100 )
		self.fileName = fileName

	def fill( self, function, zlabel, saveName ):
		for x in self.xList:
			for y in self.yList:
				if "W_gsq_" in self.fileName:
					if x == 500 and y == 520 or x == 1000 and y == 1520:
						continue # these points are not available for wino sample
				gHisto = readHisto( self.fileName, "gMet%s_%s"%(x,y) )
				eHisto = readHisto( self.fileName, "eMet%s_%s"%(x,y) )
				fHisto = readHisto( self.fileName, "fMet%s_%s"%(x,y) )
				fUpHisto = readHisto( self.fileName, "fMetUp%s_%s"%(x,y) )
				fDownHisto = readHisto( self.fileName, "fMetDown%s_%s"%(x,y) )
				if gHisto and eHisto and fHisto and fUpHisto and fDownHisto:
					eHisto = applyFakeRateEWK( eHisto )
					self.histo.SetBinContent( self.histo.FindBin(x,y), function( gHisto, eHisto, fHisto, fUpHisto, fDownHisto ) )
		self.histo.GetZaxis().SetTitle( zlabel )
		self.draw( saveName )

	def draw( self, saveName ):
		can = ROOT.TCanvas()
		can.cd()
		# When the first Palette object is created (Draw("colz")), it will inherit
		# it's title from the histogram. If the histogram is changing, the palette
		# won't recognize.
		palette = self.histo.GetListOfFunctions().FindObject("palette")
		if palette:
			palette.GetAxis().SetTitle( self.histo.GetZaxis().GetTitle() )
		self.histo.Draw("colz")
		SaveAs( can, saveName )

	def xsection( self, file_ ):
		xsec = readSignalXSection( file_ )
		for x in self.xList:
			for y in self.yList:
				self.histo.SetBinContent( self.histo.FindBin(x,y), xsec[(x,y)][0] )
		self.histo.GetZaxis().SetTitle( "#sigma [pb]" )
		self.draw( "xsectionSignal" )


def signalContamination( gHisto, eHisto, fHisto, fUpHisto, fDownHisto ):
	return 1.*addHistos([ fHisto, eHisto ]).Integral() / gHisto.Integral() if gHisto.Integral() else 0

def signalContaminationErrorRel( gHisto, eHisto, fHisto, fUpHisto, fDownHisto ):
	e1 = ROOT.Double()
	e2 = ROOT.Double()
	i1 = addHistos([ fHisto, eHisto ]).IntegralAndError(0,-1, e1 )
	i2 = gHisto.IntegralAndError(0,-1, e2 )
	return sqrt( (e1/i1)**2 + (e2/i2)**2 ) if i1 and i2 else 0

def signalContaminationError( gHisto, eHisto, fHisto, fUpHisto, fDownHisto ):
	e1 = ROOT.Double()
	e2 = ROOT.Double()
	i1 = addHistos([ fHisto, eHisto ]).IntegralAndError(0,-1, e1 )
	i2 = gHisto.IntegralAndError(0,-1, e2 )
	return i1/i2 * sqrt( (e1/i1)**2 + (e2/i2)**2 ) if i1 and i2 else 0

def signalContaminationE( gHisto, eHisto, fHisto, fUpHisto, fDownHisto ):
	return 1.*eHisto.Integral() / gHisto.Integral() if gHisto.Integral() else 0

#sp = SignalPlane( "W_gsq_01.root" )
#sp.fill( signalContamination, "signal contamination", "signalContaminationPlane" )
#sp.fill( signalContaminationError, "#sigma_{signal contermination}", "signalConmatinationError" )
#sp.fill( signalContaminationErrorRel, "#sigma_{signal cont.} / signal cont.", "signalConmatinationErrorRel" )
#sp.fill( signalContaminationE, "signal contamination e#rightarrow#gamma", "signalContaminationPlaneE" )
#sp.xsection( "/home/knut/master/infos/Spectra_gsq_W_8TeV.xsec" )



