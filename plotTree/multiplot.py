import ROOT
import ConfigParser
import os.path
from treeFunctions import *

datasetConf = ConfigParser.SafeConfigParser()
datasetConf.read( "dataset.cfg" )

integratedLumi = 19300 #pb

class Dataset:
	def __init__( self, filename, treeName, cut="1", label=None, color=None ):
		self.tree = readTree( filename, treeName )
		self.additionalCut = cut
		self.label = label

		for configName in datasetConf.sections():
			if filename.count( configName ):
				self.datasetLabel = datasetConf.get( configName, "label" )
				self.color = datasetConf.get( configName, "color" )

		if color:
			self.color = color

class Plot:
	def __init__( self ):
		pass
	# binning, etc

class SingleHisto:
	def __init__( self, name, histo, label=False ):
		self.histo = histo
		self.label = label
		self.name = name


class Multihisto:
	def __init__( self ):
		self.histos = []
		self.stack = []
		self.denominator = None
		self.numerator = None
		self.leg = myLegend(.7,.7,.95,.92)
		self.legendOption = "l"

	def addHisto( self, singleHisto, label=None, toStack=False, draw="hist" ):
		if toStack:
			self.stack.append( singleHisto )
		else:
			self.histos.append( (singleHisto, draw) )
		if label:
			if "hist" in draw:
				self.legendOption = "l"
			if draw == "e2":
				self.legendOption = "f"
			self.leg.AddEntry( singleHisto, label, self.legendOption )

	def stackHistos( self ):
		if self.stack:
			stack = ROOT.THStack()
			for h in self.stack:
				stack.Add( h )
		#self.histo.append( stack )

	def GetMinimum( self ):
		mini = 100000
		for hist,draw in self.histos:
			# search for minimum > 0
			if hist.GetMinimum(0) < mini:
				mini = hist.GetMinimum(0)
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
