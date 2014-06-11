import argparse
import ROOT

def readHisto( fileName, histoname, newName="" ):

	f = ROOT.TFile( fileName )
	ROOT.gROOT.cd()

	# This histogram is still associated with the file and will be lost after
	# the end of the function.
	histoInFile = f.Get( histoname ).Clone()
	return histoInFile


if __name__ == "__main__":

	arguments = argparse.ArgumentParser( description="Print out filters" )
	arguments.add_argument( "-o", "--outputFilename", default="puWeights.root" )
	arguments.add_argument( "--mc", default="mc.root"  )
	arguments.add_argument( "--data", default="data.root" )
	opts = arguments.parse_args()

	s7Scenario = readHisto( opts.mc, "pileupScenarioS7" )
	s10Scenario = readHisto( opts.mc, "pileupScenarioS10" )

	dataHist = readHisto( opts.data, "pileup" )
	dataHistUp = readHisto( opts.data, "pileupUp" )
	dataHistDown = readHisto( opts.data, "pileupDown" )

	for h in [ s7Scenario, s10Scenario, dataHist, dataHistUp, dataHistDown ]:
		h.Scale( 1./h.Integral() )

	pileupWeightS10 = dataHist.Clone( "pileupWeightS10" )
	pileupWeightS10Up = dataHistUp.Clone( "pileupWeightUpS10" )
	pileupWeightS10Down = dataHistDown.Clone( "pileupWeightDownS10" )

	for h in pileupWeightS10, pileupWeightS10Up, pileupWeightS10Down:
		h.Divide( s10Scenario )

	pileupWeightS7 = dataHist.Clone( "pileupWeightS7" )
	pileupWeightS7Up = dataHistUp.Clone( "pileupWeightUpS7" )
	pileupWeightS7Down = dataHistDown.Clone( "pileupWeightDownS7" )

	for h in pileupWeightS7, pileupWeightS7Up, pileupWeightS7Down:
		h.Divide( s7Scenario )

	from createMChist import writeObjectsToFile
	writeObjectsToFile( [pileupWeightS7, pileupWeightS7Up, pileupWeightS7Down, pileupWeightS10, pileupWeightS10Up, pileupWeightS10Down ], opts.outputFilename )
