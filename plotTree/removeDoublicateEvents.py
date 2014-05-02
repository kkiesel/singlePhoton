#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

def removeDublicateEvents( file, treename ):
	tree = readTree( file, treename )

	treename = "booleanTree%s"%treename
	newTree = ROOT.TTree( treename, " blub" )
	import numpy
	x1 = numpy.zeros( 1, dtype=float )
	newTree.Branch( "x", x1, "x/D" )

	json1 = []

	for e in tree:
		if (e.runNumber, e.luminosityBlockNumber, e.eventNumber) in json1:
			print "Double event in %s from %s"%( tree.GetName(), tree.GetFile().GetName() )
			x1[0] = 0
		else:
			x1[0] = 1
			json1.append( (e.runNumber, e.luminosityBlockNumber, e.eventNumber) )
		newTree.Fill()

	f = ROOT.TFile( file, "update")
	f.cd()
	newTree.Write("", ROOT.TObject.kOverwrite)
	f.Close()


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--trees", nargs="+", default=["photonTree","photonJetTree","photonElectronTree"] )
	opts = arguments.parse_args()

	for file in opts.filenames:
		for tree in opts.trees:
			removeDublicateEvents( file, tree )

