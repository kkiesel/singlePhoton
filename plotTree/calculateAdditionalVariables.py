#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from sys import stdout
from treeFunctions import *
import numpy
import math

class variableToTree:
	"""This class is used to attach variables to a tree. The variables will be saved
	and the function will be evaluated each entry of the original tree."""
	def __init__(self, tree, name, function):
		self.name = name
		self.function = function
		self.x = numpy.zeros( 1, dtype=float )
		tree.Branch( name, self.x, name+"/D" )

	def calc(self, event):
		self.x[0] = self.function( event )

def metLLFunc( e ):
	metVect = ROOT.TVector3()
	metVect.SetPtEtaPhi( e.met, 0, e.metPhi )

	for l in e.electrons:
		lVect = ROOT.TVector3()
		lVect.SetPtEtaPhi( l.pt, l.eta, l.phi )
		metVect += lVect
	for l in e.muons:
		lVect = ROOT.TVector3()
		lVect.SetPtEtaPhi( l.pt, l.eta, l.phi )
		metVect += lVect

	return metVect.Pt()

def dPhiGammaMet( e ):
	return e.photons[0].DeltaPhi( e.metPhi )

def gammaTight( photon ):
	return photon.ptJet() > 110 \
		and photon.chargedIso < 2.6 \
		and photon.neutralIso < 3.5+0.04*photon.ptJet() \
		and photon.photonIso < 1.3+0.005*photon.ptJet()

def gammaLoose( photon ):
	return photon.ptJet() > 110 \
		and ( photon.chargedIso < 5.2 or (photon.neutralIso < 3.5+0.04*photon.ptJet() and photon.photonIso < 1.3+0.005*photon.ptJet()))\
		and ( photon.neutralIso < 7+0.06*photon.ptJet() or (photon.chargedIso < 2.6 and photon.photonIso < 1.3+0.005*photon.ptJet()))\
		and ( photon.photonIso < 2.6+0.0075*photon.ptJet() or (photon.chargedIso < 2.6 and photon.neutralIso < 3.5+0.04*photon.ptJet()))\
		and not gammaTight( photon )

def recoilChristian( e ):
	recoil = ROOT.TVector3()
	tmpVect = ROOT.TVector3()

	thisPhoton = None
	if e.GetName() == "photonJetTree":
		for photon in e.photons:
			if gammaLoose( photon ):
				thisPhoton = photon
				break
	else:
		for photon in e.photons:
			if gammaTight( photon ):
				thisPhoton = photon
				break

	if not thisPhoton:
		return -10

	thisJet = e.jets.at(thisPhoton.matchedJetIndex) if thisPhoton.matchedJetIndex >=0 else thisPhoton

	for j in e.jets:
		if j.pt<30. or abs(j.eta)>3.0 or j.DeltaR( thisJet) <0.5:
			continue
		tmpVect.SetPtEtaPhi( j.pt, j.eta, j.phi )
		recoil += tmpVect
	return recoil.Pt()

def leadingGPt( e ):
	pt = -10
	if e.GetName() == "photonJetTree":
		for photon in e.photons:
			if gammaLoose( photon ):
				pt = photon.ptJet()
				break
	else:
		for photon in e.photons:
			if gammaTight( photon ):
				pt = photon.ptJet()
				break
	return pt

def M( p1, p2 ):
	# invariant mass
	e1 = ROOT.TLorentzVector()
	e2 = ROOT.TLorentzVector()
	e1.SetPtEtaPhiE( p1.pt, p1.eta, p1.phi, p1.pt/math.sin(2*math.atan(math.exp(-p1.eta))) )
	e2.SetPtEtaPhiE( p2.pt, p2.eta, p2.phi, p2.pt/math.sin(2*math.atan(math.exp(-p2.eta))) )
	return (e1+e2).M()

def mee( e ):
	if e.GetName() == "photonTree":
		if e.electrons.size() >= 2:
			return M( e.electrons.at(0), e.electrons.at(1) )
	return -10

def mmm( e ):
	if e.GetName() == "photonTree":
		if e.muons.size() >= 2:
			return M( e.muons.at(0), e.muons.at(1) )
	return -10

def mll( e ):
	if e.GetName() == "photonTree":
		if e.muons.size() >= 2:
			return M( e.muons.at(0), e.muons.at(1) )
		if e.electrons.size() >= 2:
			return M( e.electrons.at(0), e.electrons.at(1) )
	return -10

def mTMet( e, p1 ):
	e1 = ROOT.TLorentzVector()
	e2 = ROOT.TLorentzVector()
	e1.SetPtEtaPhiE( p1.pt, p1.eta, p1.phi, p1.pt/math.sin(2*math.atan(math.exp(-p1.eta))) )
	e2.SetPtEtaPhiE( e.met, 0, e.metPhi, 0 )
	return (e1+e2).Mt()

def mTe( e ):
	if e.GetName() == "photonTree":
		if e.electrons.size() >= 1:
			return mTMet( e, e.electrons.at(0) )
	return -10

def mTm( e ):
	if e.GetName() == "photonTree":
		if e.muons.size() >= 1:
			return mTMet( e, e.muons.at(0) )
	return -10

def mT( e ):
	if e.GetName() == "photonTree":
		if e.electrons.size() >= 1:
			return mTMet( e, e.electrons.at(0) )
		if e.muons.size() >= 1:
			return mTMet( e, e.muons.at(0) )
	return -10


def metZcorr( e ):
	thisMll = mll( e )
	if 60 < thisMll and thisMll < 120:
		e1 = ROOT.TLorentzVector()
		e2 = ROOT.TLorentzVector()
		met = ROOT.TLorentzVector()
		p1 = None
		p2 = None
		if e.electrons.size() >= 2:
			p1 = e.electrons.at(0)
			p2 = e.electrons.at(1)
		elif e.muons.size() >= 2:
			p1 = e.muons.at(0)
			p2 = e.muons.at(1)
		else:
			return -10
		e1.SetPtEtaPhiE( p1.pt, p1.eta, p1.phi, p1.pt/math.sin(2*math.atan(math.exp(-p1.eta))) )
		e2.SetPtEtaPhiE( p2.pt, p2.eta, p2.phi, p2.pt/math.sin(2*math.atan(math.exp(-p2.eta))) )
		met.SetPtEtaPhiE( e.met, 0, e.metPhi, e.met )
		return (met+e1+e2).Mt()
	return -10

def metWcorr( e ):
	thisMll = mT( e )
	if 50 < thisMll and thisMll < 100:
		e1 = ROOT.TLorentzVector()
		met = ROOT.TLorentzVector()
		p1 = None
		if e.electrons.size() >= 1:
			p1 = e.electrons.at(0)
		elif e.muons.size() >= 1:
			p1 = e.muons.at(0)
		else:
			return -10
		e1.SetPtEtaPhiE( p1.pt, p1.eta, p1.phi, p1.pt/math.sin(2*math.atan(math.exp(-p1.eta))) )
		met.SetPtEtaPhiE( e.met, 0, e.metPhi, e.met )
		return (met+e1).Mt()
	return -10

def e1Tagger( e ):
	#tags: (v)veto, (l)loose, (m)medium, (t)tight
	# ignore order g1g2 = g2g1
	if e.electrons.size() > 1:
		if e.electrons.at(0).isStatus( 4 ):
			return 4
		if e.electrons.at(0).isStatus( 3 ):
			return 3
		if e.electrons.at(0).isStatus( 2):
			return 2
		if e.electrons.at(0).isStatus( 1 ):
			return 1

def e2Tagger( e ):
	if e.electrons.size() > 1:
		if e.electrons.at(1).isStatus( 4 ):
			return 4
		if e.electrons.at(1).isStatus( 3 ):
			return 3
		if e.electrons.at(1).isStatus( 2):
			return 2
		if e.electrons.at(1).isStatus( 1 ):
			return 1

def createNewVariableTree( filename, treename, treeAppendix="AddVariables" ):
	import numpy

	newTreeName = treename + treeAppendix
	newTree = ROOT.TTree( newTreeName, "Tree containing additional Variables" )

	newVariables = []
	newVariables.append( variableToTree( newTree, "metLL", metLLFunc ) )
	#newVariables.append( variableToTree( newTree, "dPhiGammaMet", dPhiGammaMet ) )
	newVariables.append( variableToTree( newTree, "recoilChr", recoilChristian ) )
	newVariables.append( variableToTree( newTree, "thisPt", leadingGPt ) )

	#invariant masses
	newVariables.append( variableToTree( newTree, "mee", mee ) )
	newVariables.append( variableToTree( newTree, "mmm", mmm ) )
	newVariables.append( variableToTree( newTree, "mll", mll ) )

	newVariables.append( variableToTree( newTree, "mT", mT ) )
	newVariables.append( variableToTree( newTree, "mTe", mTe ) )
	newVariables.append( variableToTree( newTree, "mTm", mTm ) )
	newVariables.append( variableToTree( newTree, "metZcorr", metZcorr ) )
	newVariables.append( variableToTree( newTree, "metWcorr", metWcorr ) )
	newVariables.append( variableToTree( newTree, "e1Tag", e1Tagger ) )
	newVariables.append( variableToTree( newTree, "e2Tag", e2Tagger ) )

	origTree = readTree( filename, treename )
	nEvents = origTree.GetEntries()
	for event in origTree:
		#if event.eventNumber != 285841452 or event.runNumber != 190733 or event.luminosityBlockNumber != 267:
		#	continue
		if not event.GetReadEntry()%10000:
			stdout.write( "\r%s / %s"%(event.GetReadEntry(), nEvents ) )
			stdout.flush()

		for var in newVariables:
			var.calc( event )

		newTree.Fill()
	print
	return newTree


if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	opts = arguments.parse_args()

	for filename in opts.filenames:
		print filename
		f = ROOT.TFile( filename, "update" )
		f.cd()

		#for treename in [ "photonTree", "photonJetTree", "photonElectronTree" ]:
		for treename in [ "photonTree", "photonJetTree"]:
			treeFriend = createNewVariableTree( filename, treename )
			treeFriend.Write("", ROOT.TObject.kOverwrite)
		f.Close()
