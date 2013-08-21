import ROOT
from treeFunctions import myLegend

class Multihisto:
	def __init__( self ):
		self.orderByIntegral = True
		self.histos = []
		self.histosToStack = []
		self.denominator = None
		self.numerator = None
		self.leg = myLegend(.7,.7,.95,.92)
		self.legendOption = "l"

	def addHisto( self, singleHisto, label=None, toStack=False, draw="hist" ):
		if toStack:
			singleHisto.SetFillColor( singleHisto.GetLineColor() )
			singleHisto.SetLineColor( 1 )
			self.histosToStack.append( singleHisto )
			self.leg.AddEntry( singleHisto, label, "f" )
		else:
			self.histos.append( (singleHisto, draw) )
			if label:
				if "p" in draw:
					self.legendOption = "pl"
				if "hist" in draw:
					self.legendOption = "l"
				if draw == "e2":
					self.legendOption = "f"
				self.leg.AddEntry( singleHisto, label, self.legendOption )

	def stackHistos( self ):
		if not self.histosToStack:
			return
		if self.orderByIntegral:
			# sort histograms first
			self.histosToStack.sort( key=lambda x: x.Integral() )
		stack = ROOT.THStack()
		stack.SetTitle( ";%s;%s"%(self.histosToStack[0].GetXaxis().GetTitle(),self.histosToStack[0].GetYaxis().GetTitle()) )
		for h in self.histosToStack:
			stack.Add( h )
		self.histos.append( (stack,"hist") )

	def GetMinimum( self ):
		mini = 100000
		for hist,draw in self.histos:
			# search for minimum > 0
			try:
				if hist.GetMinimum(0) < mini:
					mini = hist.GetMinimum(0)
			except:
				myMin = hist.GetMinimum()
				if myMin < mini and myMin > 0:
					mini = myMin
		return mini

	def GetMaximum( self ):
		maxi = -100000
		for hist, draw in self.histos:
			if hist.GetMaximum() > maxi:
				maxi = hist.GetMaximum()
		return maxi

	def Draw( self  ):
		self.stackHistos()

		if self.histos:
			# adjust maximum and minimum for log and not log
			maximum = self.GetMaximum()
			minimum = self.GetMinimum()
			if ROOT.gPad.GetLogy():
				maximum = maximum*5
				minimum = minimum/3
			else:
				maximum = maximum + (maximum-minimum)*.1
				minimum = minimum - (maximum-minimum)*.1
			self.histos[0][0].SetMaximum( maximum )
			self.histos[0][0].SetMinimum( minimum )

			self.histos[0][0].Draw(self.histos[0][1])
			for hist, draw in self.histos[1:]:
				hist.Draw("same %s"%draw)
		else:
			print "No histogram added"

		if self.leg.GetListOfPrimitives().GetSize():
			self.leg.Draw()
