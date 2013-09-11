#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from multiplot import Multihisto
from treeFunctions import *

def drawSignalContamination( filename, xy, split ):
	import array
	label, unit, binning = readAxisConf("met")
	metBinning = array.array( "d", binning )

	gHisto = readHisto( filename, "gMet%s_%s"%xy ).Rebin(len(metBinning)-1, randomName(), metBinning )
	eHisto = readHisto( filename, "eMet%s_%s"%xy ).Rebin(len(metBinning)-1, randomName(), metBinning )
	fHisto = readHisto( filename, "fMet%s_%s"%xy ).Rebin(len(metBinning)-1, randomName(), metBinning )

	eHisto = applyFakeRateEWK( eHisto )
	eContamination = divideHistos( eHisto, gHisto )
	eContamination.SetLineColor(2)

	fContamination = divideHistos( fHisto, gHisto )
	fContamination.SetLineColor(4)

	totalContamination = divideHistos( addHistos( [eHisto, fHisto] ), gHisto )
	totalContamination.SetLineColor(1)
	totalContamination.SetMarkerColor(1)

	for h in [totalContamination, fContamination, eContamination]:
		h.GetYaxis().SetTitleOffset( 1.85 )
		h.SetTitle(";%s [%s];signal contamination"%(label,unit))
		h.SetLineWidth(2)

	mh = Multihisto()
	if split:
		mh.addHisto( totalContamination, "Total", draw="e0" )
		mh.addHisto( eContamination, "#gamma_{e}", draw="e0" )
		mh.addHisto( fContamination, "#gamma_{jet}", draw="e0" )
	else:
		mh.addHisto( totalContamination, "", draw="e0" )

	can = ROOT.TCanvas()
	can.cd()
	can.SetLogy(False)
	mh.Draw()
	SaveAs(can, "signalContamination_%s_%s_%s_%s"%(xy+(split,filename[0:-4] )))



if __name__ == "__main__":
	arguments = argparse.ArgumentParser()
	arguments.add_argument("filenames", nargs="+", type=isValidFile )
	arguments.add_argument("--xy", nargs=2, default = [1200,1220], type=int )
	arguments.add_argument("--split", action="store_true" )
	opts = arguments.parse_args()

	for inName in opts.filenames:
		drawSignalContamination( inName, tuple(opts.xy), opts.split )

