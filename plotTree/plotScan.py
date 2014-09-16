import ROOT
ROOT.gROOT.SetBatch()

def style_SUS_14_004():
	s = ROOT.TStyle( "style_SUS_14_004", "Style used for common paper" )

	s.SetPadTopMargin(.15)
	s.SetPadBottomMargin(.1)
	s.SetPadLeftMargin(.15)
	s.SetPadRightMargin(.2)

	s.SetPadTickX(1)
	s.SetPadTickY(1)

	textFont = 42
	textSize = 0.045
	s.SetLegendFont( textFont )
	s.SetTextFont( textFont )
	s.SetLabelFont( textFont )
	s.SetTextSize( textSize )
	s.SetLabelSize( textSize )

	s.SetOptStat(0)

	#s.SetOptLogz(1)
	s.SetPalette(1)
	s.SetNumberContours(999)
	s.cd()

def transparentPalette():
	from array import array
	NRGBs = 5
	NCont = 999
	stops = array("d",[0.00, 0.34, 0.61, 0.84, 1.00])
	red= array("d",[0.50, 0.50, 1.00, 1.00, 1.00])
	green = array("d",[ 0.50, 1.00, 1.00, 0.60, 0.50])
	blue = array("d",[1.00, 1.00, 0.50, 0.40, 0.50])
	ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
	ROOT.gStyle.SetNumberContours(NCont)


def readFile( filename ):
	out = {}
	f = ROOT.TFile( filename )

	ROOT.gROOT.cd()
	for prim in f.GetListOfKeys():
		name = prim.GetName()
		obj = prim.ReadObj().Clone( name+"Clone" )
		out[ name ] = obj.Clone( name )

	f.Close()
	return out

def readConfig( filename ):
	import ConfigParser
	configuration = ConfigParser.SafeConfigParser()
	configuration.read( filename )

	configDict = {}
	for sec in configuration.sections():
		configDict[sec] = {}
		for opt in configuration.options( sec ):
			configDict[sec][opt] = configuration.get( sec, opt )

	return configDict

def cut2dHisto( h, xMin=None, xMax=None, yMin=None, yMax=None ):
	xMinBin = h.GetXaxis().FindBin( xMin ) if xMin else 1
	xMaxBin = h.GetXaxis().FindBin( xMax ) if xMax else h.GetNbinsX()
	yMinBin = h.GetYaxis().FindBin( yMin ) if yMin else 1
	yMaxBin = h.GetYaxis().FindBin( yMax ) if yMax else h.GetNbinsY()

	xMin = h.GetXaxis().GetBinLowEdge( xMinBin )
	xMax = h.GetXaxis().GetBinUpEdge( xMaxBin )
	yMin = h.GetYaxis().GetBinLowEdge( yMinBin )
	yMax = h.GetYaxis().GetBinUpEdge( yMaxBin )

	newH = ROOT.TH2D( h.GetName()+"Part", h.GetTitle(), int( (xMax-xMin)/h.GetXaxis().GetBinWidth(1) ), xMin, xMax, int( (yMax-yMin)/h.GetYaxis().GetBinWidth(1) ), yMin, yMax )
	for xBin in range( h.GetNbinsX()+1 ):
		for yBin in range( h.GetNbinsY()+1 ):
			x = h.GetXaxis().GetBinCenter( xBin )
			y = h.GetYaxis().GetBinCenter( yBin )
			newH.SetBinContent( newH.FindBin(x,y), h.GetBinContent( xBin, yBin ) )

	return newH

style_SUS_14_004()
transparentPalette()

configFile = "limitPlots.cfg"

configDict = readConfig( configFile )

for name, config in configDict.iteritems():

	# Read objects from file
	f = ROOT.TFile( config["file"] )
	xsec = f.Get( config["xsec"] )

	obs = f.Get( config["obs"] )
	obsP = f.Get( config["obsp"] )
	obsM = f.Get( config["obsm"] )

	exp = f.Get( config["exp"] )
	expP = f.Get( config["expp"] )
	expM = f.Get( config["expm"] )

	if "GGM" in config["scan"]:
		xsec = cut2dHisto( xsec, xMin=500, yMin=500 )
	if "T5wg" == config["scan"]:
		xsec = cut2dHisto( xsec, xMin=800 )

	# beautify histogram
	scanObjects = ( "x", "y" )
	if "GGM" in config["scan"]:
		scanObjects = ( "#tilde{q}", "#tilde{g}" )
	elif "T5" in config["scan"]:
		charge = "0" if "gg" in config["scan"] else ""
		scanObjects = ( "#tilde{g}", "#tilde{#chi}^{%s}_{1}"%charge )
	else:
		print "Please use 'T5' or 'GGM' as scan"

	xsec.SetTitle( ";m_{%s} (GeV);m_{%s} (GeV);95%% C.L. upper limit on cross section (fb)"%scanObjects )

	xsec.SetLabelSize( ROOT.gStyle.GetTextSize(), "xyz" )
	xsec.SetTitleSize( ROOT.gStyle.GetLabelSize(), "xyz" )

	xsec.SetTitleOffset( 1.02, "x" )
	xsec.SetTitleOffset( 1.6, "y" )
	xsec.SetTitleOffset( 1.5, "z" )

	xsec.GetXaxis().SetNdivisions( 504 )
	xsec.GetYaxis().SetNdivisions( 504 )

	# set exclusion contours
	for gr in obs, obsP, obsM:
		gr.SetLineColor( 1 )
	for gr in exp, expP, expM:
		gr.SetLineColor( 2 )
		gr.SetLineStyle( 7 )
	for gr in exp, obs:
		gr.SetLineWidth( 4 )
	for gr in obsP, obsM, expP, expM:
		gr.SetLineWidth( 2 )


	can = ROOT.TCanvas( "can_"+name, "", 600, 600 )

	if "T5" in config["scan"]:
		xsec.SetMinimum(1)
		can.SetTopMargin( .06 )
		can.SetLogz()
	if "GGM" in config["scan"]:
		can.SetLogz()
		can.SetTopMargin( .15 )
	if "GGMBino" in config["scan"]:
		xsec.SetMinimum(1)
	if "GGMWino" in config["scan"]:
		xsec.SetMinimum(9.9)

	xsec.Draw("colz")

	for gr in obs, obsP, obsM, exp, expM, expP:
		gr.Draw(" same l ")


	# Draw additional information
	if "GGM" in config["scan"]:
		leg = ROOT.TLegend( .5, .86, 1, .95 )
	elif "T5" in config["scan"]:
		leg = ROOT.TLegend( .15, .84, .65, .93 )
	leg.SetBorderSize(0)
	leg.SetFillStyle(0)
	leg.AddEntry( obs, "Observed #pm 1 #sigma_{theory}", "l" )
	leg.AddEntry( exp, "Expected #pm 1 #sigma_{experiment}", "l" )
	leg.Draw()

	legShift = .01
	legUp = leg.Clone()
	legUp.Clear()
	legUp.SetY1( legUp.GetY1() + legShift )
	legUp.SetY2( legUp.GetY2() + legShift )
	legUp.AddEntry( obsP, " ", "l" )
	legUp.AddEntry( expP, " ", "l" )
	legUp.Draw()

	legDown = leg.Clone()
	legDown.Clear()
	legDown.SetY1( legDown.GetY1() - legShift )
	legDown.SetY2( legDown.GetY2() - legShift )
	legDown.AddEntry( obsM, " ", "l" )
	legDown.AddEntry( expM, " ", "l" )
	legDown.Draw()

	infoText = ""
	if config["analysis"] == "singlePhoton":
		infoText = "CMS                          19.7 fb^{-1} (8 TeV)           #geq1#gamma#geq2jets"
	elif config["analysis"] == "razor":
		infoText = "CMS                          19.7 fb^{-1} (8 TeV)           #geq2#gamma#geq1jets"
	else:
		print "analysis has to be 'singlePhoton' or 'razor'"

	info = ROOT.TLatex( 0, .96, infoText )
	info.SetNDC()
	info.Draw()

	if "GGMWino" in config["scan"]:
		model = ROOT.TLatex( 0, .89, "GGM Wino-like #tilde{#chi}^{0}_{1}, m_{#tilde{#chi}^{0}_{1}} = 375GeV" )
		model.SetTextSize( model.GetTextSize()*.8 )
	if "GGMBino" in config["scan"]:
		model = ROOT.TLatex( 0, .89, "GGM Bino-like #tilde{#chi}^{0}_{1}, m_{#tilde{#chi}^{0}_{1}} = 375GeV" )
		model.SetTextSize( model.GetTextSize()*.8 )
	elif config["scan"] == "T5gg":
		model = ROOT.TLatex( .16, .79, "pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrowqq#tilde{#chi}^{0}_{1}, #tilde{#chi}^{0}_{1}#rightarrow#gamma#tilde{G}" )
		model.SetTextSize( model.GetTextSize()*.8 )
	elif config["scan"] == "T5wg":
		model = ROOT.TLatex( .16, .79, "pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrowqq#tilde{#chi}_{1}, #tilde{#chi}_{1}#rightarrow#gamma/W#tilde{G}" )
		model.SetTextSize( model.GetTextSize()*.8 )

	model.SetNDC()
	model.Draw()

	ROOT.gPad.SaveAs( "plots/{}.pdf".format(name) )
