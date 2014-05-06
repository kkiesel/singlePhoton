#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from sys import stdout
from treeFunctions import *

def createNewVariableTree( filename, treename, treeAppendix="AddVariables" ):
	import numpy

	newTreeName = treename + treeAppendix
	newTree = ROOT.TTree( newTreeName, "Tree containing additional Variables" )

	dPhiGammaMet = numpy.zeros( 1, dtype=float )
	dPhiGammaMetxy = numpy.zeros( 1, dtype=float )
	dPhiGammaJet1 = numpy.zeros( 1, dtype=float )
	dPhiGammaJet2 = numpy.zeros( 1, dtype=float )
	dPhiMetJet1 = numpy.zeros( 1, dtype=float )
	dPhiMetJet2 = numpy.zeros( 1, dtype=float )
	dRGammaJet1 = numpy.zeros( 1, dtype=float )
	dRGammaJet2 = numpy.zeros( 1, dtype=float )
	mht = numpy.zeros( 1, dtype=float )
	recoil = numpy.zeros( 1, dtype=float )

	mhtVect = ROOT.TVector3()
	recoilVect = ROOT.TVector3()
	tmpVect = ROOT.TVector3()

	# declaration of variables
	newTree.Branch( "dPhiGammaMet", dPhiGammaMet, "dPhiGammaMet/D" )
	newTree.Branch( "dPhiGammaMetxy", dPhiGammaMetxy, "dPhiGammaMetxy/D" )
	newTree.Branch( "dPhiGammaJet1", dPhiGammaJet1, "dPhiGammaJet1/D" )
	newTree.Branch( "dPhiGammaJet2", dPhiGammaJet2, "dPhiGammaJet2/D" )
	newTree.Branch( "dPhiMetJet1", dPhiMetJet1, "dPhiMetJet1/D" )
	newTree.Branch( "dPhiMetJet2", dPhiMetJet2, "dPhiMetJet2/D" )
	newTree.Branch( "dRGammaJet1", dRGammaJet1, "dRGammaJet1/D" )
	newTree.Branch( "dRGammaJet2", dRGammaJet2, "dRGammaJet2/D" )
	newTree.Branch( "mht", mht, "mht/D" )
	newTree.Branch( "recoil", recoil, "recoil/D" )

	origTree = readTree( filename, treename )

	nEvents = origTree.GetEntries()
	for e in origTree:
		if not e.GetReadEntry()%10000:
			stdout.write( "\r%s / %s"%(e.GetReadEntry(), nEvents ) )
			stdout.flush()

		dPhiGammaMet[0] = e.photons[0].DeltaPhi( e.metPhi )
		dPhiGammaMetxy[0] = e.photons[0].DeltaPhi( e.metShiftxyPhi )

		# set default values in case not enought jets
		dPhiGammaJet1[0] = -5
		dPhiMetJet1[0] = -5
		dPhiGammaJet2[0] = -5
		dPhiMetJet2[0] = -5

		foundJet1 = False
		foundJet2 = False

		mhtVect.SetXYZ(0,0,0)
		recoilVect.SetXYZ(0,0,0)

		for j in e.jets:
			if j.isMatch( 3 ): # jet was used to cout, so good jet
				if abs(j.eta)< 2.5 and j.DeltaR( e.photons[0] ) > 0.4:
					tmpVect.SetPtEtaPhi( j.pt, j.eta, j.phi )
					mhtVect += tmpVect
					recoilVect += tmpVect

				if not foundJet1:
					dPhiGammaJet1[0] = e.photons[0].DeltaPhi( j.phi )
					dPhiMetJet1[0] = -j.DeltaPhi( e.metPhi )
					dRGammaJet1[0] = e.photons[0].DeltaR( j )
					foundJet1 = True
					continue # to look for second jet
				if foundJet1 and not foundJet2:
					dPhiGammaJet2[0] = e.photons[0].DeltaPhi( j.phi )
					dPhiMetJet2[0] = -j.DeltaPhi( e.metPhi )
					dRGammaJet2[0] = e.photons[0].DeltaR( j )
					break #, since both jets were found

		tmpVect.SetPtEtaPhi( e.photons[0].pt, e.photons[0].eta, e.photons[0].phi )
		mhtVect += tmpVect

		mht[0] = mhtVect.Pt()
		recoil[0] = recoilVect.Pt()

		newTree.Fill()
	print

	return newTree


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	opts = arguments.parse_args()

	for filename in opts.filenames:
		f = ROOT.TFile( filename, "update" )
		f.cd()

		for treename in [ "photonTree", "photonJetTree", "photonElectronTree" ]:
			treeFriend = createNewVariableTree( filename, treename )
			treeFriend.Write(ROOT.TObject.kSingleKey)
		f.Close()
