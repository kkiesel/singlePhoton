"""Here functions are collected which prettify the plots and the way of coding."""

def isValidFile( fileName ):
	"""This function can be used for argparse to check if the file exists."""
	import argparse
	import os
	if not os.path.exists( fileName ):
		raise argparse.ArgumentTypeError( "File %s does not exist"%fileName )
	else:
		return fileName

def getDatasetAbbr( fileName ):
	"""Parses the abbrevation for a sample from a root fileName."""
	import re
	match = re.match(".*slim(.*)_V.*_tree.root", fileName )
	if match:
		return match.groups()[0]
	else:
		return fileName

def datasetToLatex( datasetAbbr ):
	"""Translates the dataset name to a TLatex name"""
	sets = { "AllQCD": "(#gamma+)QCD",
			"GJets": "#gamma+QCD",
			"TTbar": "t#bar{t}",
			"WJet": "W"
			}
	for part, label in sets.iteritems():
		if part in datasetAbbr:
			return label
	return datasetAbbr

def createDatasetLabel( datasetAbbr ):
	import ROOT
	"""Creates sample info which then can be printed."""
	datasetLabel = ROOT.TPaveText(.4,.94,.6,1, "ndc")
	datasetLabel.SetFillColor(0)
	datasetLabel.SetBorderSize(0)
	datasetLabel.AddText( datasetToLatex(datasetAbbr) )
	return datasetLabel

def randomName():
	"""
	Generate a random string. This function is useful to give ROOT objects
	different names to avoid overwriting.
	"""
	from random import randint
	from sys import maxint
	return "%x"%(randint(0, maxint))

def roundToSignificantDigits(x, sig=2):
	"""Round number to 'sig' significant digits. If the number is large enough,
	just print the integer.
	"""
	from math import log10, floor
	if x >= 10**(sig-1):
		return int(round(x))
	if x>0:
		return round(x, sig-int(floor(log10(x)))-1)
	elif x<0:
		return round(-x, sig-int(floor(log10(-x)))-1)
	return x

def myLegend( x1, y1, x2=0,y2=0 ):
	import ROOT
	if x2==0 or y2==0:
		style = ROOT.gROOT.GetStyle("tdrStyle")
		x2 = 1 - style.GetPadRightMargin()
		y2 = 1 - style.GetPadTopMargin()
	leg = ROOT.TLegend(x1,y1,x2,y2)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)
	return leg

def readAxisConf( plot, configurationFileName="axis.cfg" ):
	"""Read the configuration file for the axis.
	returns the label, the unit and the binning as list if avaible
	"""
	import ConfigParser
	configuration = ConfigParser.SafeConfigParser()
	configuration.read( configurationFileName )
	#brackets are identified as sections, so they have to be deleted
	plot = plot.replace("[","").replace("]","")
	if not configuration.has_section( plot ):
		return "","",""
	label = configuration.get( plot, "label" )
	unit = configuration.get( plot, "unit" )
	binning = configuration.get( plot, "binning" )
	if binning:
		binning = map(float, binning.split(" "))
	else:
		binning = []
	return label, unit, binning

def manipulateSaveName( saveName ):
	"""Replace some charakters, so root nor unix have problems to read them."""
	#saveName = saveName.replace("/","VS")
	saveName = saveName.replace(" ","_")
	unallowedCharacters = ["{","}","(",")","#","|",".","[","]","/","$"]
	for char in unallowedCharacters:
		saveName = saveName.replace( char, "" )
	return saveName

def SaveAs( can, folder, name, endings=["pdf"] ):
	"""Save ROOT.TCanvas in specified folder with a cleaned plot name."""
	for ending in endings:
		can.SaveAs( folder+"/"+manipulateSaveName( name )+"."+ending )

class PlotCaption:
	"""Creates the superscription for each plot, eg
	'19fb^{-1} sqrt{s)=8TeV #geq1#gamma #geq2jets'
	"""
	def __init__( self, x0=.96, y0=.96, analysisInfo=True, option="ndc" ):
		import ROOT
		self.x0 = x0
		self.text = ROOT.TLatex( x0, y0, "" )
		self.text.SetTextSize(0.03)
		self.text.SetNDC()
		if analysisInfo:
			self.addAnalysisInfo()

	def addAnalysisInfo( self, lumi=19800, e=8, defaultcuts="#geq1#gamma,#geq2jets" ):
		self.text.SetText( self.text.GetX(), self.text.GetY(), "%.1ffb^{-1} #sqrt{s}=%sTeV %s"%( lumi/1000., e, defaultcuts ) )

	def appendEnd( self, string ):
		newText = self.text.GetTitle() + string
		self.text.Clear()
		self.text.SetText( self.text.GetX(), self.text.GetY(), newText )

	def appendFront( self, string ):
		newText = string + self.text.GetTitle()
		self.text.Clear()
		self.text.SetText( self.text.GetX(), self.text.GetY(), newText )

	def controlCut( self ):
		self.appendEnd(",#slash{E}_{T}<100GeV")

	def signalCut( self ):
		self.appendEnd(",#slash{E}_{T}#geq100GeV")

	def Draw( self ):
		import ROOT
		shiftNDC = self.text.GetXsize() / ( ROOT.gPad.GetX2() - ROOT.gPad.GetX1() )
		self.text.SetX( self.x0-shiftNDC )
		self.text.Draw()


