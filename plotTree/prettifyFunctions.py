"""Here functions are collected which prettify the plots and the way of coding."""

def isValidFile( fileName ):
	"""This function can be used for argparse to check if the file exists."""
	import argparse
	import os
	if not os.path.isfile( fileName ):
		raise argparse.ArgumentTypeError( "File %s does not exist"%fileName )
	else:
		return fileName

def getDatasetAbbr( fileName, slim=True ):
	"""Parses the abbrevation for a sample from a root fileName."""
	import re

	datamatch = re.match("PhotonHad([ABCD]{0,1})_V.*", fileName.split("/")[-1])
	if datamatch:
		return "Data %s"%datamatch.groups()[0]

	prefix = "slim" if slim else ""
	match = re.match("%s(.*)_V.*_tree.root"%prefix, fileName.split("/")[-1] )
	if match:
		return match.groups()[0]
	else:
		return fileName

def replaceListElementsBy( list_, search, replacement ):
	if len([x for x in list_ if x in search ]) == len(search):
		return [x for x in list_ if x not in search ]+[replacement]
	return list_

def mergeDatasetAbbr( datasetAbbrs ):
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["GJets_200_400", "GJets_400_inf"], "GJets" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["QCD_250_500", "QCD_500_1000", "QCD_1000_inf"], "QCD" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["GJets", "QCD"], "AllQCD" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["TTJets", "WJets"], "EWK" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["WGamma", "ZGamma"], "FSR" )
	return datasetAbbrs


def datasetToLatex( datasetAbbr ):
	"""Translates the dataset name to a TLatex name"""
	sets = { "AllQCD": "(#gamma+)QCD",
			"GJets": "#gamma+QCD",
			"TTJets": "t#bar{t}",
			"WJets": "W",
			"QCD": "QCD"
			}
	#for part, label in sets.iteritems():
	#	if part in datasetAbbr:
	#		return label
	try:
		return sets[datasetAbbr]
	except:
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

def getLumiWeight( datasetAbbr, nGenEvents, integratedLumi=19800, configName="dataset.cfg" ):
	import ConfigParser
	datasetConf = ConfigParser.SafeConfigParser()
	datasetConf.read( configName )
	if datasetConf.has_option( datasetAbbr, "crosssection" ):
		crosssection = datasetConf.getfloat( datasetAbbr, "crosssection" )
	else:
		raise NameError( "Configuration for %s not found"%datasetAbbr )

	return 1. * integratedLumi * crosssection / nGenEvents

def manipulateSaveName( saveName ):
	"""Replace some charakters, so root nor unix have problems to read them."""
	#saveName = saveName.replace("/","VS")
	saveName = saveName.replace(" ","_")
	unallowedCharacters = ["{","}","(",")","#","|",".","[","]","/","$"]
	for char in unallowedCharacters:
		saveName = saveName.replace( char, "" )
	return saveName

def SaveAs( can, name, folder="plots", endings=["pdf"] ):
	"""Save ROOT.TCanvas in specified folder with a cleaned plot name."""
	for ending in endings:
		can.SaveAs( folder+"/"+manipulateSaveName( name )+"."+ending )

class PlotCaption:
	"""Creates the superscription for each plot, eg
	'19fb^{-1} sqrt{s)=8TeV #geq1#gamma #geq2jets'
	"""
	def __init__( self, x0=.96, y0=.96, analysisInfo=True, option="ndc", signal=False, control=False ):
		import ROOT
		self.x0 = x0
		self.text = ROOT.TLatex( x0, y0, "" )
		self.text.SetTextSize(0.03)
		self.text.SetNDC()
		if analysisInfo:
			self.addAnalysisInfo()
		if signal:
			self.signalCut()
		if control:
			self.controlCut()

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

def readSignalXSection( filename ):
	"""Read xsection and other informations for various signal MC from a file
	found at https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA3PrivateSignalMC2012
	The syntax is printf("%6.0f %6.0f %6.0f %6.0f %6.0f LO: %9.3e + %9.3e - %9.3e NLO: %9.3e + %9.3e - %9.3e\n", $nevents, $msquark, $mgluino, $mbino, $mwwino, $loxsec,$lohierr,$loloerr,$nloxsec,$nlohierr,$nloloerr);

	returns list[ (point1, point2) ] = ( sigmaNLO, errorUp, errorDown )
	"""
	f = open( filename )
	text = f.readlines()
	f.close()

	info = {}

	floatMatch = "\-{0,1}\d\.\d+e[+-]\d{2}"
	import re
	for t in text:
		matches = re.match(" (\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+) LO: (%s) \+ (%s) \- (%s) NLO: (%s) \+ (%s) \- (%s)\n"%(floatMatch,floatMatch,floatMatch,floatMatch,floatMatch,floatMatch), t ).groups()

		info[ (int(matches[1]), int(matches[2]) ) ] = ( float(matches[8]), float(matches[9]), float(matches[10]) )

	return info

