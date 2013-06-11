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
	if not isinstance(nBins, int):
		result.Scale(1,"width")
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

def appendOverflowBin( oldHist, overflow ):
	import ROOT
	import array
	oldArray = oldHist.GetXaxis().GetXbins()
	newArray = array.array('d', [0]*(oldArray.GetSize()+1))
	for i in range( oldArray.GetSize() ):
		newArray[i] = oldArray[i]
	newArray[-1] = newArray[-2] + overflow
	newHist = ROOT.TH1F( randomName(), "", len(newArray)-1, newArray )

	for bin in range(1, len(newArray)):
		newHist.SetBinContent(bin, oldHist.GetBinContent(bin ))
		newHist.SetBinError(bin, oldHist.GetBinError(bin ))
	return newHist

def extractHisto( dataset, plot, overflow=0 ):
	label, unit, binning = readAxisConf( plot )
	histo = createHistoFromTree( dataset.tree, plot, "weight*(%s)"%(dataset.additionalCut), nBins=binning)
	if overflow > 0:
		histo = appendOverflowBin(histo, overflow)

	histo.SetLineColor( dataset.color )
	histo.SetMarkerColor( dataset.color )

	ytitle = "Entries"
	if binning:
		ytitle+= " / Bin"
		if unit:
			label+= " [%s]"%unit
	else:
		ytitle+= " / %s"%histo.GetBinWidth(1)
		if unit:
			label+= " [%s]"%unit
			ytitle+= " %s"%unit
	histo.SetTitle(";%s;%s"%(label, ytitle))
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

def readAxisConf( plot, configurationFileName="axis.cfg" ):
	import ConfigParser
	configuration = ConfigParser.SafeConfigParser()
	configuration.read( configurationFileName )
	#brackets are identified as sections, so they have to be deleted
	plot = plot.replace("[","").replace("]","")
	label = configuration.get( plot, "label" )
	unit = configuration.get( plot, "unit" )
	binning = configuration.get( plot, "binning" )
	binning = map(int, binning.split(" "))
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
