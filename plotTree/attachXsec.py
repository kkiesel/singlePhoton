#! /usr/bin/env python2
# -*- coding: utf-8 -*-
from sys import stdout
from treeFunctions import *
import re


def modify( inputFileName, printOnly ):
	"""Key function of splitCanidates.py. The main loop and object selection is
	defined here."""
	datasetAbbr = getDatasetAbbr( inputFileName, slim=False )
	if "Data" in datasetAbbr:
		return

	if inputFileName.startswith("slim"):
		print "File %s already processed?"%inputFileName
		return

	eventHisto = readHisto( inputFileName )
	processNEvents = eventHisto.GetBinContent(1)
	lumiWeight = getLumiWeight( datasetAbbr, processNEvents )
	if printOnly:
		return datasetAbbr, processNEvents, lumiWeight

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
		if not inTree.GetEntries():
			continue
		tree = inTree.CloneTree(0)
		ROOT.gROOT.ProcessLine("struct variablesToChange { Float_t weight(1); };")
		change = ROOT.variablesToChange()
		tree.SetBranchAddress("weight", ROOT.AddressOf( change, "weight" ) )

		for event in inTree:
			change.weight = lumiWeight * event.weight
			tree.Fill()
		fout.cd()
		tree.AutoSave("overwrite")

	for histName in histNames:
		h = readHisto( inputFileName, histName )
		h.Scale( lumiWeight )
		fout.cd()
		h.Write()

	fout.Close()

	return datasetAbbr, processNEvents, lumiWeight

if __name__ == "__main__":

	arguments = argparse.ArgumentParser( description="Slim tree" )
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--printOnly", action="store_true" )
	opts = arguments.parse_args()

	rawTable = []
	for inName in opts.filenames:
		out = modify( inName, opts.printOnly )
		if out:
			abbr, nGen, w = out
			s = w*nGen/19712
			rawTable.append( (abbr, s, nGen ) )

	print
	print "Sample name  &  $\sigma$ [pb]  &  nGen [1e6]"
	for abbr, s, nGen in rawTable:
		niceAbbr = abbr.replace("GJets","#gammaJets").replace("TTJets","t#bar{t}").replace("ZGamma", "#gammaZ")
		print "\\verb|%s|  &  %.3f  &  %s \\\\"%(abbr, s, nGen/1e6 )

