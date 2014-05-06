#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *

def removeDublicateEvents( filename, treename, printOutput ):
	oldTree = readTree( filename, treename )

	newTree = oldTree.CloneTree(0)

	json1 = []

	for e in oldTree:
		if (e.runNumber, e.luminosityBlockNumber, e.eventNumber) not in json1:
			newTree.Fill()
			json1.append( (e.runNumber, e.luminosityBlockNumber, e.eventNumber) )
		elif printOutput:
			print "Double event %s,%s,%s in %s from %s"%( e.runNumber, e.luminosityBlockNumber, e.eventNumber, oldTree.GetName(), oldTree.GetFile().GetName() )

	f = ROOT.TFile( filename, "update")
	f.cd()
	newTree.Write("", ROOT.TObject.kOverwrite)
	f.Close()


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--trees", nargs="+", default=["photonTree","photonJetTree","photonElectronTree"] )
	arguments.add_argument( "-v", "--verbose", action="store_true" )
	opts = arguments.parse_args()

	for filename in opts.filenames:
		for tree in opts.trees:
			removeDublicateEvents( filename, tree, opts.verbose )

