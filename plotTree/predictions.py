from treeFunctions import *

def multiDimFakeRate( filenames, plot="met", cut="1", isData=True ):

	if isData:
		# from AN2013_340_v2
		weightString = "(1. - 0.993 * (1. - std::pow(photons[0].pt / 2.9 + 1., -2.4)) * (1. - 0.23 * std::exp(-0.2777 * nTracksPV))* (1. - 5.66e-4 * nVertex))"
		# relative uncertainty: 11%
	else:
		# yutaro's email, 01.05.2014
		weightString = "1. - 0.998 * (1. - std::pow(photons[0].pt / 8.67 + 1., -3.93)) * (1. - 0.198 * std::exp(-0.4 * nTracksPV)) * (1. - 1.27e-4 * nVertex)"

	eHist = None
	for filename in filenames:
		eTree = readTree( filename, "photonElectronTree" )
		#try:
		#	eTree.AddFriend( "photonElectronTreeangular", filename )
		#except:
		#	pass

		recE = getHisto( eTree, plot, weight="weight*(%s)"%weightString, color=2, fillEmptyBins=True, cut=cut )
		recE.SetFillColor( recE.GetLineColor() )

		if eHist:
			eHist.Add( recE )
		else:
			eHist = recE

	#recE = addRelativeUncertainty( recE, 2.75892/100 )
	return eHist


def getMixedWeigthHisto( filenames, predFilenames, commonCut, control=True ):
	"""Calculate #photons/#photonFakes in bins of photons.ptJet and a second
	(global) variable.

	filenames: files containing photons
	predFilenames: files containing fakes
	"""

	regionCut = "met<100" if control else "met>=100"

	xVar = "photons[0].ptJet()"
	yVar = "ht"
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

	numerator = None
	for fileName in filenames:
		gTree = readTree( fileName, "photonTree" )
		num = createHistoFromTree2D( gTree, yVar+":"+xVar, "weight*( %s && %s )"%(regionCut, commonCut), xbinning, ybinning )
		if numerator:
			numerator.Add( num )
		else:
			numerator = num

	denominator = None
	for fileName in predFilenames:
		foTree = readTree( fileName, "photonJetTree" )
		den = createHistoFromTree2D( foTree, yVar+":"+xVar, "weight*( %s && %s )"%(regionCut, commonCut), xbinning, ybinning )
		if denominator:
			denominator.Add( den )
		else:
			denominator = den

	weight2D = divideHistos( numerator, denominator )

	# Set the weight and error for empty bins above the diagonal to one.
	if False:
		for i in range( weight2D.GetXaxis().GetNbins()+1 ):
			for j in range( weight2D.GetYaxis().GetNbins()+1 ):
				if not weight2D.GetBinContent( i, j ):
					weight2D.SetBinContent( i, j, 1 )
					weight2D.SetBinError( i, j, 0 )

	return weight2D


def getTreeFriendFromWeights( tree, h_weight, weightTreeName ):
	"""Write weight for a tree in another tree in a given file.
	This tree can be added to the original tree via 'AddFriend()'.

	fileName: name of file to which tree is written
	tree: tree which is weighted
	h_weight: two dimensional histogram with weights
	weighTreeName: name of the new tree
	"""

	weightTree = ROOT.TTree( weightTreeName, "Tree containing QCD weights" )
	import numpy
	weight = numpy.zeros( 1, dtype=float)
	weight_error = numpy.zeros( 1, dtype=float)
	weightTree.Branch( "w_qcd", weight, "w_qcd/D" )
	weightTree.Branch( "w_qcd_error", weight_error, "w_qcd_error/D" )

	from sys import stdout
	for event in tree:
		if not event.GetReadEntry()%10000:
			stdout.write( "\r%s / %s"%(event.GetReadEntry(), event.GetEntries() ) )
			stdout.flush()

		b = h_weight.FindBin( event.photons.at(0).ptJet(), event.ht )
		weight[0] = h_weight.GetBinContent( b )
		weight_error[0] = h_weight.GetBinError( b )
		weightTree.Fill()
	print

	return weightTree

def attachWeightsToFiles( filenames, weight2D, weightTreeName ):
	for fileName in filenames:
		foTree = readTree( fileName, "photonJetTree" )
		treeFriend = getTreeFriendFromWeights( foTree, weight2D, weightTreeName )

		f = ROOT.TFile( fileName, "update" )
		f.cd()
		treeFriend.Write( "", ROOT.TObject.kSingleKey )
		f.Close()


def predictionHistos( filenames, plot, cut, modifyEmptyBins ):
	fHist, sysHist = None, None

	for filename in filenames:
		fTree = readTree( filename, "photonJetTree" )
		fTree.AddFriend( "foWeights", filename )

		# Prediction + statistical uncertainty coming from the loose control region
		hist = getHisto( fTree, plot, weight="weight*w_qcd", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )

		# This histogram contains the uncertainty due to the weight for small met
		# Here, the bin content is the uncertainty, the bin error is meaningless.
		sHist = getHisto( fTree, plot, weight="weight*w_qcd_error", cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )

		# You can also do an approximation to calculate the error of the weight:
		doWeightApproximation = False
		if doWeightApproximation:
			sHistUp = getHisto( fTree, plot, weight="weight*(w_qcd_error+w_qcd)",
					cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
			sHistDown = getHisto( fTree, plot, weight="weight*(w_qcd-w_qcd_error)",
					cut=cut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
			for bin in range( sHistUp.GetNbinsX()+2 ):
				sHistUp.SetBinContent( bin, 0.5*(sHistUp.GetBinContent(bin)-sHistDown.GetBinContent(bin)))
			sHist = sHistUp
		# Approximation end ###################################################

		if fHist:
			fHist.Add( hist )
			sysHist.Add( sHist )
		else:
			fHist = hist
			sysHist = sHist

	fHist.SetMarkerSize(0)
	sysHist.SetFillColor( sysHist.GetLineColor() )
	sysHist.SetLineColor( sysHist.GetLineColor() )
	sysHist.SetFillStyle(3254)
	sysHist.SetMarkerSize(0)

	return fHist, sysHist

