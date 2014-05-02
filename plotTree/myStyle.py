import ROOT

def texStyle():
	style = ROOT.TStyle("texStyle", "Style to generate PGF/TikZ")

	# Canvas
	style.SetCanvasColor( ROOT.kWhite )
	style.SetCanvasBorderMode(0)
	style.SetCanvasDefH(1000)
	style.SetCanvasDefW(2000)

	# Pad:
	style.SetPadBorderMode(0)

	# Margins:
	style.SetPadTopMargin(0.05)
	style.SetPadBottomMargin(0.13)
	style.SetPadLeftMargin(0.16)
	style.SetPadRightMargin(0.02)

	# axis titles
	style.SetTitleFont(42, "xyz")
	style.SetTitleSize(0.09, "xyz")
	style.SetLabelSize(0.09, "xyz")

	# Statistic box
	style.SetOptStat(0)

	# Other
	style.SetPalette(1)
	# x-value is approximation of text width, y value large
	style.SetPaperSize(14.65,50.)

	style.SetPadTickX(1)
	style.SetPadTickY(1)






	style.cd()
	return style
