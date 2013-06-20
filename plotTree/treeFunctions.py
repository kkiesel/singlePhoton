def randomName():
	"""
	Generate a random string. This function is useful to give ROOT objects
	different names to avoid overwriting.
	"""
	from random import randint
	from sys import maxint
	return "%x"%(randint(0, maxint))

def createHistoFromTree(tree, variable, weight="", nBins=20, firstBin=None, lastBin=None ):
	"""
	tree: tree to create histo from
	variable: variable to plot (must be a branch of the tree)
	weight: weights to apply (e.g. "var1*(var2 > 15)" will use weights from var1 and cut on var2 > 15
	nBins, firstBin, lastBin: number of bins, first bin and last bin (same as in TH1F constructor)
	nBins: if nBins is a list, and to a int, a user binned plot will be generated
	returns: histogram
	"""
	from ROOT import TH1F
	name = randomName()
	if isinstance( nBins, list ):
		import array
		xBins = array.array('d', nBins )
		result = TH1F(name, variable, len(nBins)-1, xBins)
		result.Sumw2()
		tree.Draw("%s>>%s"%(variable, name), weight, "goff")
		result.Scale(1,"width")
	elif firstBin==None and lastBin==None:
		import ROOT
		tree.Draw("%s>>%s"%(variable,name), weight, "goff")
		result = ROOT.gDirectory.Get( name )
		result.Sumw2() # applying the errors here is perhaps not entirely correct
	else:
		result = TH1F(name, variable, nBins, firstBin, lastBin)
		result.Sumw2()
		tree.Draw("%s>>%s"%(variable, name), weight, "goff")
	return result

def readTree( filename, treename = "susyTree" ):
	"""
	filename: name of file containing the tree
	treename: name of the tree
	returns: TChain Object
	"""
	import ROOT
	tree = ROOT.TChain( treename )
	tree.AddFile( filename )
	return tree

def readHisto( filename, histoname="eventNumbers" ):
	"""
	filename: name of file containing the histogram
	histoname: name of the histogram
	returns: Object with the given name, should be a histogram
	"""
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
	"""Append the overflow bin to a histogram.
	oldHist: input histogram
	overflow: size of the overflow bin, which should be drawn
	"""
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

def getXMinXMax( histo_list ):
	"""Find minimum and maximum on the x-axis for a list of histograms.
	histo_lists: list of histograms
	"""
	from sys import maxint
	mini = maxint
	maxi = -maxint
	for histo in histo_list:
		mini = min( mini, histo.GetBinLowEdge(1) )
		lastBin = histo.GetNbinsX()
		maxi = max( maxi, histo.GetBinLowEdge(lastBin)+histo.GetBinWidth(lastBin) )
	return mini, maxi

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

def getHisto( tree, plot, cut="1", overflow=0, weight="weight", color=1, nBins=20, firstBin=None, lastBin=None ):
	"""Creates a histogram and apply the axis settings
	cut: cutstring applied to the tree
	overflow: size of the overflow bin, if overflow>0
	"""
	label, unit, binning = readAxisConf( plot )
	if binning:
		histo = createHistoFromTree( tree, plot, "%s*(%s)"%(weight, cut), nBins=binning)
	else:
		histo = createHistoFromTree( tree, plot, "%s*(%s)"%(weight, cut), nBins=nBins, firstBin=firstBin, lastBin=lastBin )
	if overflow > 0:
		histo = appendOverflowBin(histo, overflow)

	histo.SetLineColor( color )
	histo.SetMarkerColor( color )
	histo.SetLineWidth(2)

	ytitle = "Entries"
	if not label:
		import re
		objVarExpr = "([a-zA-Z]+)\[{0,1}(\d*)\]{0,1}\.([a-zA-Z]+)" # matches eg photon.pt
		if "Length$(" in plot:
			obj, nObj, var = re.match( "Length\$\(%s\)"%objVarExpr, plot ).groups()
			label = "N_{%s}"%obj

		elif "." in plot:
			obj, nObj, var = re.match( objVarExpr, plot ).groups()
			var = var.replace("phi","#phi_{")
			var = var.replace("eta","#eta_{")
			var = var.replace("pt","p_{T ")
			obj = obj.replace("gamma","#gamma")
			obj = obj.replace("electron","e")
			obj = obj.replace("muon", "#mu")
			if nObj:
				label = "%s%s.%s}"%(var,int(nObj)+1,obj)
			else:
				label = "%s%s}"%(var,obj)
		else:
			label = plot
	if binning:
		ytitle+= " / Bin"
		if unit:
			label+= " [%s]"%unit
	else:
		if histo.GetBinWidth(1) != 1:
			ytitle+= " / {}".format(roundToSignificantDigits(histo.GetBinWidth(1),2))
		if unit:
			label+= " [%s]"%unit
			ytitle+= " %s"%unit
	histo.SetTitle(";%s;%s"%(label, ytitle))
	return histo

def getQCDErrorHisto( tree, plot, cut="1", overflow=0, nBins=20, firstBin=None, lastBin=None ):
	"""Applies w_qcd +- w_qcd_error for qcd error propagation.
	The returned histo's content is the mean, the error are the shifts up and down.
	"""
	label, unit, binning = readAxisConf( plot )
	if binning:
		histoUp = createHistoFromTree( tree, plot, "(weight*(w_qcd+w_qcd_error))*(%s)"%(cut), nBins=binning)
		histoDown = createHistoFromTree( tree, plot, "(weight*(w_qcd-w_qcd_error))*(%s)"%(cut), nBins=binning)
	else:
		histoUp = createHistoFromTree( tree, plot, "(weight*(w_qcd+w_qcd_error))*(%s)"%(cut), nBins=nBins, firstBin=firstBin, lastBin=lastBin )
		histoDown = createHistoFromTree( tree, plot, "(weight*(w_qcd-w_qcd_error))*(%s)"%(cut), nBins=nBins, firstBin=firstBin, lastBin=lastBin )
	if overflow > 0:
		histoUp = appendOverflowBin(histoUp, overflow)
		histoDown = appendOverflowBin(histoDown, overflow)

	outHisto = histoUp.Clone( randomName() )
	for bin in range( histoUp.GetNbinsX()+1 ):
		up = histoUp.GetBinContent(bin)
		down = histoDown.GetBinContent(bin)
		outHisto.SetBinContent( bin, (up+down)/2 )
		outHisto.SetBinError( bin, (up-down)/2 )
	return outHisto

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
	"""Divide numerator by denominator."""
	option = "B" if bayes else ""
	resultHisto = numerator.Clone( randomName() )
	resultHisto.Divide( numerator, denominator, 1,1, option )
	return resultHisto

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
