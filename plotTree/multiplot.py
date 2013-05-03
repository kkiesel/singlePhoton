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
		self.drawRatio = True
		self.leg = myLegend(.7,.7,.95,.92)

	def addHisto( self, singleHisto, label=None, toStack=False ):
		if toStack:
			self.stack.append( singleHisto )
		else:
			self.histos.append( singleHisto )
		self.leg.AddEntry( singleHisto, label, "l" )

	def stackHistos( self ):
		if self.stack:
			stack = ROOT.THStack()
			for h in self.stack:
				stack.Add( h )
		#self.histo.append( stack )

	def GetMinimum( self ):
		mini = 100000
		for hist in self.histos:
			if hist.GetMinimum() < mini:
				mini = hist.GetMinimum(0)
		if ROOT.gROOT.GetStyle("tdrStyle").GetOptLogy():
			return 1
		else:
			return mini

	def GetMaximum( self ):
		maxi = -100000
		for hist in self.histos:
			if hist.GetMaximum() > maxi:
				maxi = hist.GetMaximum()
		if ROOT.gROOT.GetStyle("tdrStyle").GetOptLogy():
			return maxi*5
		else:
			return maxi*1.2

	def Draw( self, option ):
		self.stackHistos()

		if self.histos:
			self.histos[0].SetMaximum( self.GetMaximum() )
			self.histos[0].SetMinimum( self.GetMinimum() )

			self.histos[0].Draw( option )
			for hist in self.histos[1:]:
				hist.Draw("same %s"%option)
		else:
			print "No histogram added"

		if self.leg.GetListOfPrimitives().GetSize():
			self.leg.Draw()
