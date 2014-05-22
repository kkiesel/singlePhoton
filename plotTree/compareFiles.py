#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *


def getEventNumberDiff( tree1, tree2 ):
	json1 = []
	json2 = []

	for e in tree1:
		if (e.runNumber, e.luminosityBlockNumber, e.eventNumber) in json1:
			print "double event"
		break
		json1.append( (e.runNumber, e.luminosityBlockNumber, e.eventNumber) )

	for e in tree2:
		if (e.runNumber, e.luminosityBlockNumber, e.eventNumber) in json2:
			print "double event"
		json2.append( (e.runNumber, e.luminosityBlockNumber, e.eventNumber) )

	res1 = [ x for x in json1 if x not in json2 ]
	res2 = [ x for x in json2 if x not in json1 ]

	return res1, res2

def getTreeFriendsWithBooleanVariable( tree, diffList ):

	treename = "booleanTree%s"%len(diffList)
	newTree = ROOT.TTree( treename, " blub" )
	import numpy
	x1 = numpy.zeros( 1, dtype=float )
	newTree.Branch( "x", x1, "x/D" )

	for e in tree:
		x1[0] = 1 if (e.runNumber, e.luminosityBlockNumber, e.eventNumber) in diffList else 0
		newTree.Fill()

	filename = "test%s.root"%len(diffList)
	f = ROOT.TFile( filename, "recreate")
	f.cd()
	newTree.Write()
	f.Close()
	tree.AddFriend( treename, filename  )

	return tree

def getVersion( filename ):
	import re
	m = re.match( ".*(V\d\d\.\d\d)_tree.root", filename )
	if m:
		return m.groups()[0]
	else:
		return filename.replace("_tree.root", "")

def compareFiles( plot, filenames, treename, diff ):
	tree1 = readTree( filenames[0], treename )
	tree2 = readTree( filenames[1], treename )

	if diff:
		cut1, cut2 = getEventNumberDiff( tree1, tree2 )
		tree1 = getTreeFriendsWithBooleanVariable( tree1, cut1)
		tree2 = getTreeFriendsWithBooleanVariable( tree2, cut2)
		h1 = createHistoFromTree( tree1, plot, "x" )
		h2 = createHistoFromTree( tree2, plot, "x" )

	else:

		h1 = getHisto( tree1, plot, weight="1" )
		h2 = getHisto( tree2, plot, weight="1" )
		mini = min(h1.GetBinLowEdge(1), h2.GetBinLowEdge(1))
		maxi = max(h2.GetBinLowEdge(h2.GetNbinsX()+2), h1.GetBinLowEdge(h1.GetNbinsX()+2))
		h1 = getHisto( tree1, plot, weight="1", firstBin=mini, lastBin=maxi  )
		h2 = getHisto( tree2, plot, weight="1", firstBin=mini, lastBin=maxi )


	h2.SetLineColor(2)
	mh = Multihisto()
	name1 = getVersion(filenames[0])
	name2 = getVersion(filenames[1])
	mh.addHisto( h1, name1, draw="hist e" )
	mh.addHisto( h2, name2, draw="hist e" )
	mh.Draw()
	from myRatio import Ratio
	r = Ratio( "%s/%s"%(name1,name2), h1, h2 )
	r.draw()

	plot ="test"
	ROOT.gPad.GetCanvas().SaveAs( "compareFiles_%s.pdf"%plot )




if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs=2, type=isValidFile )
	arguments.add_argument( "--plots", nargs="+", default=[ "met"] )
	arguments.add_argument( "--diff", action="store_true" )
	arguments.add_argument( "--tree", default="photonTree" )
	opts = arguments.parse_args()

	for plot in opts.plots:
		compareFiles( plot, opts.filenames, opts.tree, opts.diff )

