from treeFunctions import *

fCut = " && (Max$(photons.chargedIso) < 5.2 || (Max$(photons.neutralIso-0.04*photons.pt) < 3.5 && Max$(photons.photonIso-0.005*photons.pt) < 1.3)) \
&& (Max$(photons.neutralIso-0.06*photons.pt) < 7 || (Max$(photons.chargedIso) < 2.6 && Max$(photons.photonIso-0.005*photons.pt) < 1.3 )) \
&& (Max$(photons.photonIso-0.0075*photons.pt) < 2.6 || (Max$(photons.chargedIso) < 2.6 && Max$(photons.neutralIso-0.04*photons.pt) < 3.5))"

#fCut = " && (Max$(photons.chargedIso) < 5.2 || (Max$(photons.neutralIso-0.04*photons.ptJet()) < 3.5 && Max$(photons.photonIso-0.005*photons.ptJet()) < 1.3)) \
#&& (Max$(photons.neutralIso-0.06*photons.ptJet()) < 7 || (Max$(photons.chargedIso) < 2.6 && Max$(photons.photonIso-0.005*photons.ptJet()) < 1.3 )) \
#&& (Max$(photons.photonIso-0.0075*photons.ptJet()) < 2.6 || (Max$(photons.chargedIso) < 2.6 && Max$(photons.neutralIso-0.04*photons.ptJet()) < 3.5))"

#fCut = " && Max$( (photons[0].chargedIso/2.6)**2 + ((photons[0].neutralIso-0.04*photons[0].pt)/3.5)**2 + ((photons[0].photonIso-0.005*photons[0].pt)/1.3)**2 < 20 )" # eleptic cut
#fCut = " && 1"

#gCut = " && Max$(photons.chargedIso)<0.7 && Max$(photons.neutralIso-0.04*photons.pt)< 0.4 && Max$(photons.photonIso-0.005*photons.pt) < 0.5 & Max$(photons.sigmaIetaIeta<)0.011" #tight working point
gCut = " && photons.chargedIso<2.6 && photons.neutralIso-0.04*photons.ptJet()< 3.5 && photons.photonIso-0.005*photons.ptJet() < 1.3" # loose working point with pt*
#gCut = "&&1"
#fCut = "&&1"
yVar = "recoil"


def multiDimFakeRate( filenames, plot="met", cut="1", isData=True ):

	if isData:
		# from AN2013_340_v2
		weightString = "(1. - 0.993 * (1. - std::pow(photons[0].pt / 2.9 + 1., -2.4)) * (1. - 0.23 * std::exp(-0.2777 * nTracksPV))* (1. - 5.66e-4 * nVertex))"
		# relative uncertainty: 11%
	else:
		# yutaro's email, 01.05.2014, this formula is crab
		weightString = "1. - 0.998 * (1. - std::pow(photons[0].pt / 8.67 + 1., -3.93)) * (1. - 0.198 * std::exp(-0.4 * nTracksPV)) * (1. - 1.27e-4 * nVertex)"
		# yutaro's email, 12.02.2014
		weightString = "1. - 0.995003 * (1. - TMath::Power(1.968e-01 * photons[0].pt + 1., -3.120e+00)) * (1. - 4.392e-01 * TMath::Exp(-3.394e-01 * nTracksPV)) * (1. - 2.308e-04 * nVertex)"
		# yutaro's email, 06.05.2014
		weightString = "1 - (1 - 0.00623) * (1 - std::pow(photons[0].pt / 4.2 + 1,-2.9)) * (1 - 0.29 * std::exp(-0.335 * nTracksPV)) * (1 - 0.000223 * nVertex)"


	eHist = None
	for filename in filenames:
		eTree = readTree( filename, "photonElectronTree" )
		recE = getHisto( eTree, plot, weight="weight*(%s)"%weightString, color=2, fillEmptyBins=True, cut=cut )

		if eHist:
			eHist.Add( recE )
		else:
			eHist = recE

	return eHist


def getMixedWeigthHisto( filenames, predFilenames, commonCut, control=True, fillEmptyBins=False ):
	"""Calculate #photons/#photonFakes in bins of photons.ptJet and a second
	(global) variable.

	filenames: files containing photons
	predFilenames: files containing fakes
	"""

	regionCut = "met<100" if control else "met>=100"

	xVar = "photons[0].ptJet()"
	xlabel, xunit, xbinning = readAxisConf( xVar )
	ylabel, yunit, ybinning = readAxisConf( yVar )

	numerator = None
	for fileName in filenames:
		gTree = readTree( fileName, "photonTree" )
		num = createHistoFromTree2D( gTree, yVar+":"+xVar, "weight*( %s && %s )"%(regionCut, commonCut+gCut), xbinning, ybinning )
		if numerator:
			numerator.Add( num )
		else:
			numerator = num

	denominator = None
	for fileName in predFilenames:
		foTree = readTree( fileName, "photonJetTree" )
		den = createHistoFromTree2D( foTree, yVar+":"+xVar, "weight*( %s && %s )"%(regionCut, commonCut+fCut), xbinning, ybinning )
		if denominator:
			denominator.Add( den )
		else:
			denominator = den

	weight2D = divideHistos( numerator, denominator )

	weightIntegral = numerator.Integral() / denominator.Integral()

	# Set the weight to the global weight
	if fillEmptyBins:
		for i in range( weight2D.GetXaxis().GetNbins()+1 ):
			for j in range( weight2D.GetYaxis().GetNbins()+1 ):
				if not weight2D.GetBinContent( i, j ):
					weight2D.SetBinContent( i, j, weightIntegral )
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

		b = h_weight.FindBin( event.photons.at(0).ptJet(), eval("event.%s"%yVar) )
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
		treeFriend.Write( "", ROOT.TObject.kOverwrite )
		f.Close()


def predictionHistos( filenames, plot, cut, modifyEmptyBins=True, modifyEmptyWeightBins=False ):
	fHistSum, sysHistSum, sysHistEmptyBinSum = None, None, None

	for filename in filenames:
		fTree = readTree( filename, "photonJetTree" )
		fTree.AddFriend( "foWeights", filename )

		# Prediction + statistical uncertainty coming from the loose control region
		fHist = getHisto( fTree, plot, weight="weight*w_qcd", cut=cut+fCut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )

		# This histogram contains the uncertainty due to the weight for small met
		# Here, the bin content is the uncertainty, the bin error is meaningless.
		sysHist = getHisto( fTree, plot, weight="weight*w_qcd_error", cut=cut+fCut, color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )

		if modifyEmptyWeightBins:
			sysHistEmptyBin = getHisto( fTree, plot, weight="weight", cut="w_qcd*(%s)+(w_qcd<0.0001)*(%s)"%(cut,cut), color=46, firstBin=1, lastBin=1, fillEmptyBins=modifyEmptyBins )
			sysHistEmptyBin.Add( hist, -1. ) # get difference to normal prediction
		else:
			sysHistEmptyBin = None

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

		if fHistSum:
			fHistSum.Add( fHist )
			sysHistSum.Add( sysHist )
		else:
			fHistSum = fHist
			sysHistSum = sysHist
		if sysHistEmptyBin:
			sysHistEmptyBinSum.Add( sysHistEmptyBin )
		else:
			sysHistEmptyBinSum = sysHistEmptyBin

	fHist.SetMarkerSize(0)
	sysHist.SetFillColor( sysHist.GetLineColor() )
	sysHist.SetLineColor( sysHist.GetLineColor() )
	sysHist.SetFillStyle(3254)
	sysHist.SetMarkerSize(0)

	for bin in range(sysHistSum.GetNbinsX()+2):
		uncert = sysHistSum.GetBinContent(bin) |qPlus| sysHistEmptyBinSum.GetBinContent(bin) if modifyEmptyWeightBins else sysHistSum.GetBinContent(bin)
		sysHistSum.SetBinError( bin, uncert )
		sysHistSum.SetBinContent( bin, fHistSum.GetBinContent(bin) )

	return fHistSum, sysHistSum

