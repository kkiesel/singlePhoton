#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *


def getEventNumberDiff( tree1, tree2 ):
	json1 = []
	json2 = []

	for e in tree1:
		if (e.runNumber, e.luminosityBlockNumber, e.eventNumber) in json1:
			print "Double event in %s from %s"%( tree1.GetName(), tree1.GetFile().GetName() )
		json1.append( (e.runNumber, e.luminosityBlockNumber, e.eventNumber) )

	for e in tree2:
		if (e.runNumber, e.luminosityBlockNumber, e.eventNumber) in json2:
			print "Double event in %s from %s"%( tree2.GetName(), tree2.GetFile().GetName() )
		json2.append( (e.runNumber, e.luminosityBlockNumber, e.eventNumber) )

	res1 = [ x for x in json1 if x not in json2 ]
	res2 = [ x for x in json2 if x not in json1 ]

	print "Events in t1:", len(json1)
	print "only in t1:  ", len(res1)

	print "Events in t2:", len(json2)
	print "only in t2:  ", len(res2)

	return res1, res2



if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument( "filenames", nargs="+", type=isValidFile )
	arguments.add_argument( "--treenames", nargs="+", default=["photonTree"] )
	opts = arguments.parse_args()

	nFiles = len(opts.filenames)
	nTrees = len(opts.treenames)

	if nFiles == 1:
		files = opts.filenames*2
	elif nFiles == 2:
		files = opts.filenames
	elif nFiles > 2:
		raise NameError("To many filenames provided")
	else:
		raise NameError("files ???")

	if nTrees == 1:
		treenames = opts.treenames*2
	elif nTrees == 2:
		treenames = opts.treenames
	elif nFiles > 2:
		raise NameError("To many treenames provided")
	else:
		raise NameError("trees ???")

	t1 = readTree( files[0], treenames[0] )
	t2 = readTree( files[1], treenames[1] )

	r1, r2 = getEventNumberDiff( t1, t2 )


