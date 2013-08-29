#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from sys import stdout
from treeFunctions import *
import re

def modify( inputFileName, processNEvents, genMatch=False ):
	if genMatch:
		print "not implemented yet"
	"""Key function of splitCanidates.py. The main loop and object selection is
	defined here."""
	datasetAbbr = getDatasetAbbr( inputFileName, slim=False )
	print "Processing file %s with %s configuration"%(inputFileName, datasetAbbr)

	eventHisto = readHisto( inputFileName )
	if processNEvents < 0:
		processNEvents = eventHisto.GetBinContent(1)
	lumiWeight = getLumiWeight( datasetAbbr, processNEvents )

	treeNames = []
	histNames = []

	f = ROOT.TFile( inputFileName )
	for item in f.GetListOfKeys():
		if item.GetName() == "eventNumbers":
			continue # This histogram is used to scale
		if item.GetClassName() == "TTree":
			treeNames.append( item.GetName() )
		elif re.match("TH[123].", item.GetClassName() ):
			histNames.append( item.GetName() )
		else:
			print "Object %s could not be orderer in one of the classes."%item.GetName()

	# Output tree definition
	import os
	outputFileName = "slim"+os.path.basename( inputFileName )
	fout = ROOT.TFile( outputFileName, "recreate" )

	for treeName in treeNames:
		inTree = readTree( inputFileName, treeName )
		tree = inTree.CloneTree(0)
		ROOT.gROOT.ProcessLine("struct variablesToChange { Float_t weight(1); };")
		change = ROOT.variablesToChange()
		tree.SetBranchAddress("weight", ROOT.AddressOf( change, "weight" ) )

		for event in inTree:
			change.weight = lumiWeight * event.weight
			tree.Fill()
		fout.cd()
		tree.Write()

	for histName in histNames:
		h = readHisto( inputFileName, histName )
		h.Scale( lumiWeight )
		fout.cd()
		h.Write()

	fout.Close()

if __name__ == "__main__":

	arguments = argparse.ArgumentParser( description="Slim tree" )
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--test", action="store_true" )
	arguments.add_argument("--genMatch", action="store_true" )
	opts = arguments.parse_args()

	# set limit for number of events for testing reason
	processNEvents = 10000 if opts.test else -1

	for inName in opts.filenames:
		modify( inName, processNEvents, opts.genMatch )

