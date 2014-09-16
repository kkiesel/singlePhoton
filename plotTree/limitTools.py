import ROOT
import re
from math import sqrt

from prettifyFunctions import randomName
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

def removeUnnecessaryWhitespaces( string ):
	# remove leading, trailing and multibles of whitespaces
	string = re.sub( "^\s+", "", string ) # remove leading
	string = re.sub( "\s+$", "", string ) # remove trailing
	# replace groups of whitespaces with a single space
	string = re.sub( "\s+", " ", string )
	return string

def removeComments( string ):
	# remove all characters after a "#"
	return re.sub( "#.*", "", string )

def getFilesFromFolder( path, regex="*" ):
	from glob import glob
	return glob( "{}/{}".format( path, regex ) )

def fileToDict( filename ):
	f = open( filename )
	lines = f.readlines()
	f.close()
	content = {}
	for line in lines:
		line = removeUnnecessaryWhitespaces( line )
		line = removeComments( line )
		group = re.split( "\s*=\s*", line )
		if len(group) == 2:
			first = group[0]
			if group[1].isdigit():
				second = int( group[1] )
			else:
				try: second = float( group[1] )
				except:
					try: second = [ float(i) for i in group[1].split( " " ) if i ]
					except: second = group[1]
			if second:
				content[ first ] = second
	return content

def interpolateEmptyBins( histo, minNeighbours=3 ):
	"""A point is filled, if minNeighbours points are filled."""
	fillings = {}
	for x in range( 1, histo.GetNbinsX()+1 ):
		for y in range( 1, histo.GetNbinsY()+1 ):
			if histo.GetBinContent(x,y) == 0:
				values = [ histo.GetBinContent( x, y+1 ),
						histo.GetBinContent( x, y-1 ),
						histo.GetBinContent( x+1, y ),
						histo.GetBinContent( x-1, y ) ]
				values = [ v for v in values if v > 0 ]
				if len(values) >= minNeighbours:
					fillings[(x,y)] = sum(values)/len(values)

	for (x,y), value in fillings.iteritems():
		histo.SetBinContent( x, y, value )

	return histo

def getHistoFromScan( pointList ):
	xList = sorted(set(zip(*pointList)[0]))
	yList = sorted(set(zip(*pointList)[1]))
	stepSizeX = xList[1] - xList[0]
	stepSizeY = yList[1] - yList[0]

	if xList != range( xList[0], xList[-1]+1, stepSizeX ):
		print "Error: no equidistant x-binning"
		print xList
		print "Assume T5wg sample and set bins manually"
		xList = range( xList[0], xList[-1]+1, 50 )
	if yList != range( yList[0], yList[-1]+1, stepSizeY ):
		print "Error: no equidistant y-binning"
		print yList

	return ROOT.TH2F( randomName(), "", len(xList), xList[0]-0.5*stepSizeX, xList[-1]+0.5*stepSizeX, \
		len(yList), yList[0]-0.5*stepSizeY, yList[-1]+0.5*stepSizeY )

def getPositionOfPoint( x, y, relWidth=.2 ):
	"""A 2d histogram can be divided into 1 middle rectangle and 4 trapezoids:
	center, left, right, top, bottom.
	"""

	hist = [ i for i in ROOT.gPad.GetListOfPrimitives() if isinstance( i, ROOT.TH2 ) or isinstance( i, ROOT.TGraph2D ) ][0]

	xMin = hist.GetXaxis().GetXmin()
	xMax = hist.GetXaxis().GetXmax()
	yMin = hist.GetYaxis().GetXmin()
	yMax = hist.GetYaxis().GetXmax()

	if x > xMax or x < xMin or y > yMax or y < yMin:
		print "Warning: Point not in axis range!"
		return False

	dxMin = ( x-xMin ) / ( xMax-xMin )
	dxMax = ( xMax-x ) / ( xMax-xMin )
	dyMin = ( y-yMin ) / ( yMax-yMin )
	dyMax = ( yMax-y ) / ( yMax-yMin )

	minDist = min( dxMin, dyMin, dxMax, dyMax )

	if minDist > relWidth:
		return "center"
	if minDist == dxMin:
		return "left"
	if minDist == dxMax:
		return "right"
	if minDist == dyMin:
		return "bottom"
	if minDist == dyMax:
		return "top"

def smoothGraph( graph, n=3, sigmaCorrection=1 ):
	"""The points of TGraph 'graph' are shifted slightly, to smooth the curve.
	The next 'n' points are used from both directions, and are considered with
	a gaussian. The width of the gaussian can be set by 'sigmaCorrection'.
	Points close to the boarder are treated seperatly.
	"""

	newGraph = graph.Clone( randomName() )
	nGraph = newGraph.GetN()

	if n > nGraph: n = nGraph

	sigma = n*sigmaCorrection
	lim = 3*n
	fb = ROOT.TF1("fb", "gaus(0)", -lim, lim)
	fb.SetParameter( 0, 1 )
	fb.SetParameter( 1, 0 )
	fb.SetParameter( 2, sigma )

	gaus = [ fb.Eval(i) for i in range(-n, n+1) ]
	gaus = [ i/sum(gaus) for i in gaus ]

	for iPoint in range( nGraph ):
		# point iPoint will be changed
		x  = ROOT.Double(0)
		y  = ROOT.Double(0)
		newGraph.GetPoint( iPoint, x, y )
		x0, y0 = float(x), float(y)

		newX = 0
		newY = 0

		if iPoint-n < 0 or iPoint+n > nGraph-1:
			# normal mode, use +- n points
			for i, iTestPoint in enumerate(range( iPoint-n, iPoint+n+1 )):
				graph.GetPoint( iTestPoint, x, y )
				newX += x * gaus[i]
				newY += y * gaus[i]
		else:
			pointPosition = getPositionOfPoint( x, y )
			if pointPosition == "center":
				# use up to n points, but the same number up and down
				modN = min( n, iPoint, nGraph-iPoint )
				for i, iTestPoint in enumerate(range( iPoint-n, iPoint+n+1 )):
					graph.GetPoint( iTestPoint, x, y )
					newX += x * gaus[i]
					newY += y * gaus[i]
			else:
				for i, iTestPoint in enumerate(range( iPoint-n, iPoint+n+1 )):
					graph.GetPoint( iTestPoint, x, y )
					newX += x * gaus[i]
					newY += y * gaus[i]
				if pointPosition in [ "left", "right" ]:
					newY = y0
				else:
					newX = x0

		newGraph.SetPoint( iPoint, x, y )
	return newGraph


def fillScanHisto( scan, name ):
	h = getHistoFromScan( scan.keys() )
	h.SetName( name )
	for point, dict in scan.iteritems():
		bin = h.FindBin( *point )
		if name in dict:
			h.SetBinContent( bin, dict[name] )
	h = interpolateEmptyBins( h )

	return h

def setPointAtBeginning( cont, x, y ):
	newCont = cont.Clone( randomName() )
	newCont.SetPoint( 0, x, y )
	x = ROOT.Double(0)
	y = ROOT.Double(0)
	for i in range( cont.GetN() ):
		cont.GetPoint( i, x, y )
		newCont.SetPoint( i+1, x, y )

	return newCont


def extrapolateToEdge( cont, hist ):
	"""Extrapolates the first and last point of a graph to the boundaries of a 2d histogram, if the points are near the edge
	Points already on the edge are skipped."""

	nGraph = cont.GetN()
	x = ROOT.Double(0)
	y = ROOT.Double(0)

	for lastPoint, next2lastPoint in [ (0, 1), ( nGraph, nGraph-1 ) ]:
		cont.GetPoint( lastPoint, x, y )
		pointPosition = getPositionOfPoint( x, y )
		x0, y0 = float(x), float(y)
		cont.GetPoint( next2lastPoint, x, y )
		x1, y1 = float(x), float(y)

		if pointPosition == "bottom":
			newY = hist.GetYaxis().GetXmin()
			if y0 in [ y1, newY ]: continue
			newX = ( newY*(x1-x0)-y0*x1+y1*x0 ) / ( y1-y0 )
		elif pointPosition == "top":
			newY = hist.GetYaxis().GetXmax()
			if y0 in [ y1, newY ]: continue
			newX = ( newY*(x1-x0)-y0*x1+y1*x0 ) / ( y1-y0 )
		elif pointPosition == "left":
			newX = hist.GetXaxis().GetXmin()
			if x0 in [ x1, newX ]: continue
			newY = ( y1*(newX-x0) - y0*(newX-x1) ) / ( x1 - x0 )
		elif pointPosition == "right":
			newX = hist.GetXaxis().GetXmax()
			if x0 in [ x1, newX ]: continue
			newY = ( y1*(newX-x0) - y0*(newX-x1) ) / ( x1 - x0 )

		if lastPoint == 0:
			cont = setPointAtBeginning( cont, newX, newY )
		else:
			cont.SetPoint( nGraph, newX, newY )

	return cont


def getGraph( scan, name ):
	h = fillScanHisto( scan, name )
	gr2D = ROOT.TGraph2D( h )
	gr2D.Draw("colz")
	ROOT.gPad.Update()
	cont = gr2D.GetContourList(1)[0]

	#cont = smoothGraph( cont )
	#cont = extrapolateToEdge( cont, h )
	return cont


def extractHistos( path, regex="*.txt.result.txt" ):
	name = path.split('/')[-1][17:]
	scan = {}
	for filename in getFilesFromFolder( path, regex ):
		thisDict = fileToDict( filename )
		if "CLs observed" not in thisDict:
			continue
		if "T5" in name:
			scanStrX = "gluino"
			scanStrY = "chi1"
		if "GMSB_SqGl" in name:
			scanStrX = "squark"
			scanStrY = "gluino"

		# Calculate variables
		# ObsR = CLs observed, etc
		try:
			theoUncert = sqrt( thisDict["signal.scale.uncertainty"]**2 + thisDict["signal.PDF.uncertainty"]**2 )
		except:
			theoUncert = sqrt( 0.05**2 + thisDict["signal.PDF.uncertainty"]**2 )
		thisDict["ObsXsecLimit"] = thisDict["Xsection.NLO"] * thisDict["CLs observed"] *1000 # convert to fb
		thisDict["ObsR"] = thisDict["CLs observed"]
		thisDict["ObsRTheoUp"] = thisDict["CLs observed"] * ( 1+theoUncert )
		thisDict["ObsRTheoDo"] = thisDict["CLs observed"] * ( 1-theoUncert )

		scan[ (thisDict[scanStrX], thisDict[scanStrY] ) ] = thisDict


	xsec = fillScanHisto( scan, "ObsXsecLimit" )


	obs = getGraph(scan, "ObsR")
	obsM = getGraph(scan, "ObsRTheoUp")
	obsP = getGraph(scan, "ObsRTheoDo")

	exp = getGraph(scan, "CLs expected")
	expM = getGraph(scan, "CLs expected m1sigma")
	expP = getGraph(scan, "CLs expected p1sigma")

	obs.SetName( "graph_Obs" )
	obsM.SetName( "graph_ObsM" )
	obsP.SetName( "graph_ObsP" )
	exp.SetName( "graph_Exp" )
	expM.SetName( "graph_ExpM" )
	expP.SetName( "graph_ExpP" )

	draw = False
	if draw:
		xsec.Draw("colz")
		for i in obs, obsM, obsP, exp, expP, expM:
			i.Draw("same * l")

		ROOT.gPad.SetLogz(1)

		ROOT.gPad.SaveAs("plots/raw_%s.pdf"%name)

	out = ROOT.TFile( "%s.root"%name, "recreate" )
	for i in xsec, obs, obsM, obsP, exp, expP, expM:
		i.Write()
	out.Close()

limitDir = "../../limits/"
extractHistos( limitDir+"2014-08-06-12-39-GMSB_SqGl_met-Wino" )
extractHistos( limitDir+"2014-08-06-12-34-GMSB_SqGl_met-Bino" )
extractHistos( limitDir+"2014-08-06-13-34-SMS_T5gg" )
extractHistos( limitDir+"2014-08-06-13-34-SMS_T5wg" )
