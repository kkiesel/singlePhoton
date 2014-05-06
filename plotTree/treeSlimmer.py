#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from sys import stdout
from treeFunctions import *





if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	opts = arguments.parse_args()

	for filename in opts.filenames:

		f = ROOT.TFile( filename, "update" )
		f.cd()

		tree =  readTree( filename, "photonJetTree" )
		t = tree.CopyTree( "Max$(photons.sigmaIetaIeta)<0.012 && Max$(photons.chargedIso)<26 && Min$(photons.chargedIso) > 0.26 && Max$(photons.neutralIso-0.4*photons.pt) < 35 && Min$(photons.neutralIso-0.004*photons.pt) > 0.35 && Max$(photons.photonIso-0.05*photons.pt) < 13 && Min$(photons.photonIso-0.0005*photons.pt) > 0.13" )
		t.Write("", ROOT.TObject.kOverwrite )
		f.Close()



