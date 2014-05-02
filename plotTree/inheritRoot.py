import ROOT
import math

class Legend( ROOT.TLegend ):

	def __init__( self, x1=0, y1=0, x2=0, y2=0, header="", option="brNCD" ):

		if x2==0 and y2==0:
			style = ROOT.gROOT.GetStyle("tdrStyle")
			x2 = 1 - style.GetPadRightMargin()
			y2 = 1 - style.GetPadTopMargin()


		ROOT.TLegend.__init__( self )
		self.SetFillColor(0)
		self.SetBorderSize(0)

class Bin:
	def __init__( self, histogram, bin ):
		self.content = histogram.GetBinContent(bin)
		self.error = histogram.GetBinError(bin)

class H1F( ROOT.TH1F ):
	#TODO: test all functions
	def MergeOverflow():
		nBins = self.GetNbinsX()
		self.SetBinContent( nBins, self.GetBinContent(nBins)+self.GetBinContent(nBins+1) )
		self.fSumw2.fArray[nBins] += self.fSumw2.fArray[nBins+1]
		# TODO: set bin content 0


	def __iter__(self):
		self.currentBin = 0
		#TODO: define bin, and all functions like getBinContent
		# including under and overflow
		#return range( self.GetNbinsX()+2 )
		return self

	def next(self):
		if self.currentBin > self.GetNbinsX()+1:
			raise StopIteration
		else:
			self.currentBin += 1
			return Bin( self, self.currentBin-1 )

	def GetContent(self):
		return self.GetBinContent(self.currentBin)

	#todo: systematic uncertainty

class Tree( ROOT.TChain ):
	# addfriend


if __name__ == "__main__":
	h = H1F("nide", "dinel", 10, 0, 1)
	h.FillRandom("gaus", 100)

	for bin in h:
		print bin.content
	print

	for bin in range( h.GetNbinsX()+2):
		print h.GetBinContent(bin),
