import ROOT

def readHisto( fileName, histoname, newName ):

	f = ROOT.TFile( fileName )
	ROOT.gROOT.cd()

	# This histogram is still associated with the file and will be lost after
	# the end of the function.
	histoInFile = f.Get( histoname )

	# Change the name (only temporaly), since for TH1::Clone a different name is expected
	histoInFile.SetName(histoname+"Clone")

	# Now copy the histogram with its' original name
	histo = histoInFile.Clone( newName )

	return histo


if __name__ == "__main__":

	mcFileName = "PU_dist_S10.root"
	dataFileName = "PU_dist.root"
	dataFileNameUp = "PU_dist_up.root"
	dataFileNameDown = "PU_dist_down.root"

	histName = "pileup"


	mcHist = readHisto( mcFileName, histName, "mcHisto" )
	dataHist = readHisto( dataFileName, histName, "weight" )
	dataHistUp = readHisto( dataFileNameUp, histName, "weightUp" )
	dataHistDown = readHisto( dataFileNameDown, histName, "weightDown" )

	for h in [ mcHist, dataHist, dataHistUp, dataHistDown ]:
		h.Scale( 1./h.Integral() )

	dataHist.Divide( mcHist )
	dataHistUp.Divide( mcHist )
	dataHistDown.Divide( mcHist )

	weightFile = ROOT.TFile("puWeights.root", "recreate")
	weightFile.cd()

	dataHist.Write()
	dataHistUp.Write()
	dataHistDown.Write()

	weightFile.Close()

