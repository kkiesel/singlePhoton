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

	if firstBin == None:
		firstBin = tree.GetMinimum( variable )
	if lastBin == None:
		lastBin = tree.GetMaximum( variable )

	#make a random name you could give something meaningfull here,
	#but that would make this less readable
	name = "%x"%(randint(0, maxint))
	if isinstance(nBins, int):
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
	plot = plot.replace("[","").replace("]","")
	try:
		label = configuration.get( plot, "label" )
		unit = configuration.get( plot, "unit" )
		if unit != "":
			unit = " ["+unit+"]"
	except:
		label = ""
		unit = ""
		print "Please specify %s in your axis configuration file."%plot
	return label, unit
