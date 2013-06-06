def randomName():
	from random import randint
	from sys import maxint
	return "%x"%(randint(0, maxint))


def createHistoFromTree(tree, variable, weight="", nBins=100, firstBin=None, lastBin=None, nEvents=-1):
	"""
	tree: tree to create histo from
	variable: variable to plot (must be a branch of the tree)
	weight: weights to apply (e.g. "var1*(var2 > 15)" will use weights from var1 and cut on var2 > 15
	nBins, firstBin, lastBin: number of bins, first bin and last bin (same as in TH1F constructor)
	nEvents: number of events to process (-1 = all)
	returns: histogram
	"""
	from ROOT import TH1F
	from random import randint
	from sys import maxint
	if nEvents < 0:
		nEvents = maxint

	#make a random name you could give something meaningfull here,
	#but that would make this less readable
	name = "%x"%(randint(0, maxint))
	if isinstance(nBins, int):
		if firstBin == None:
			firstBin = tree.GetMinimum( variable )
		if lastBin == None:
			lastBin = tree.GetMaximum( variable )
		result = TH1F(name, variable, nBins, firstBin, lastBin)
	else:
		# assume nBins is list
		import array
		xBins = array.array('d', nBins )
		result = TH1F(name, variable, len(nBins)-1, xBins)

	result.Sumw2()
	tree.Draw("%s>>%s"%(variable, name), weight, "goff", nEvents)
	return result

def readTree( filename, treename = "susyTree" ):
	"""
	filename: name of file containing the tree
	treename: name of the tree
	returns: TTree Object
	"""
	#import os
	#if not os.path.isfile(filename):
	#	print( "File %s does not exist"%filename)
	import ROOT
	tree = ROOT.TChain( treename )
	tree.AddFile( filename )

	return tree

def readHisto( filename, histoname="eventNumbers" ):
	import ROOT
	import os
	if not os.path.isfile(filename):
		print( "File %s does not exist"%filename)
	f = ROOT.TFile( filename )

	# with this command, the histo is cloned to root and not deleted after end
	# of function
	ROOT.gROOT.cd()
	histo = f.Get( histoname )

	# the name +Clone is only temporaly, since for TH1::Clone a different name is expected
	histo.SetName(histoname+"Clone")
	histo = histo.Clone( histoname )
	return histo

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

def readAxisConf( plot, configuration ):
	#brackets are identified as sections, so they have to be deleted
	plot = plot.replace("[","").replace("]","")
	try:
		label = configuration.get( plot, "label" )
	except:
		label = ""
		print "Please specify %s in your axis configuration file."%plot
	try:
		unit = configuration.get( plot, "unit" )
	except:
		unit = ""
	try:
		binning = configuration.get( plot, "binning" )
		binning = map(int, binning.split(" "))
	except:
		binning = []
	return label, unit, binning

def addHistos( histos, scales=None ):
	"""
	add several histos with different scales to one single histo
	histos: list of histograms
	scales: list of scales (of same size as histos)

	return: single histogram
	"""
	if not scales:
		scales = [1]*len(histos)
	else:
		if len(histos) != len(scales):
			print "Histos and scales have to have same dimension"

	sumHist = histos[0].Clone( randomName() )
	sumHist.Scale( scales[0] )

	for i, h in enumerate(histos[1:]):
		sumHist.Add( h, scales[i+1] )

	return sumHist

def divideHistos( numerator, denominator, bayes=False ):
	option = ""
	if bayes:
		option = "B"
	resultHisto = numerator.Clone( randomName() )
	resultHisto.Divide( numerator, denominator, 1,1, option )
	return resultHisto
