import ROOT
from math import sqrt
from prettifyFunctions import randomName

class Ratio:
	def __init__( self, title, numerator, denominator, sysHisto=None ):
		self.title = title
		self.numerator = numerator
		self.denominator = denominator
		self.sysHisto = sysHisto
		self.ratio = numerator.Clone( randomName() )
		self.ratioSys = denominator.Clone( randomName() )

	def draw( self, yMin=None, yMax=None ):
		for bin in range(self.numerator.GetNbinsX()+2):
			if not self.denominator.GetBinContent(bin):
				continue
			self.ratio.SetBinContent( bin, self.numerator.GetBinContent(bin) / self.denominator.GetBinContent(bin) )
			self.ratio.SetBinError( bin, self.numerator.GetBinError(bin) / self.denominator.GetBinContent(bin) )
			if self.sysHisto:
				self.ratioSys.SetBinContent( bin, self.sysHisto.GetBinContent(bin) / self.denominator.GetBinContent(bin) )
				combinedError = sqrt( self.sysHisto.GetBinError(bin)**2 + self.denominator.GetBinError(bin)**2 )
				self.ratioSys.SetBinError( bin, combinedError / self.denominator.GetBinContent(bin) )
			else:
				self.ratioSys.SetBinContent( bin, 1 )
				self.ratioSys.SetBinError( bin, self.denominator.GetBinError(bin) / self.denominator.GetBinContent(bin) )

		if yMin == None:
			yMin = min( max(0,self.ratio.GetMinimum()), .5 )
		if yMax == None:
			yMax = min( max(1.5, self.ratio.GetMaximum()), 50 )

		self.ratio.GetYaxis().SetNdivisions(2, 0, 2)
		self.ratio.SetTitleOffset(0.3, "Y")
		self.ratio.SetTitleSize(0.15, "Y")
		self.ratio.SetYTitle( self.title )
		self.ratio.SetXTitle( "" )
		self.ratio.GetXaxis().SetLabelSize(0.0)
		self.ratio.GetYaxis().SetLabelSize(0.1)
		self.ratio.SetMinimum( yMin )
		self.ratio.SetMaximum( yMax )

		self.ratioSys.SetFillStyle(3254)
		self.ratioSys.SetMarkerSize(0)
		self.ratioSys.SetFillColor(self.ratioSys.GetLineColor())

		oneLine = ROOT.TLine( self.ratio.GetBinLowEdge(1), 1.0, self.ratio.GetBinLowEdge(self.ratio.GetNbinsX()+1), 1.0)
		oneLine.SetLineStyle(2)

		return self.ratio, self.ratioSys, oneLine
