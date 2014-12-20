import ROOT
from ROOT import TStyle

def tdrStyle():
	tdrStyle = TStyle("tdrStyle","Style for P-TDR")

	# For the canvas:
	tdrStyle.SetCanvasBorderMode(0)
	tdrStyle.SetCanvasColor(ROOT.kWhite)
	tdrStyle.SetCanvasDefH(800)  #Height of canvas
	tdrStyle.SetCanvasDefW(800)  #Width of canvas
	tdrStyle.SetCanvasDefX(0)	#POsition on screen
	tdrStyle.SetCanvasDefY(0)

	# For the Pad:
	tdrStyle.SetPadBorderMode(0)
	#tdrStyle.SetPadBorderSize(Width_t size = 1)
	tdrStyle.SetPadColor(ROOT.kWhite)
	tdrStyle.SetPadGridX(False)
	tdrStyle.SetPadGridY(False)
	tdrStyle.SetGridColor(0)
	tdrStyle.SetGridStyle(3)
	tdrStyle.SetGridWidth(1)

	# For the frame:
	tdrStyle.SetFrameBorderMode(0)
	tdrStyle.SetFrameBorderSize(1)
	tdrStyle.SetFrameFillColor(0)
	tdrStyle.SetFrameFillStyle(0)
	tdrStyle.SetFrameLineColor(1)
	tdrStyle.SetFrameLineStyle(1)
	tdrStyle.SetFrameLineWidth(1)

	# For the histo:
	#tdrStyle.SetHistFillColor(1)
	#tdrStyle.SetHistFillStyle(0)
	tdrStyle.SetHistLineColor(1)
	tdrStyle.SetHistLineStyle(0)
	tdrStyle.SetHistLineWidth(2)
	#tdrStyle.SetLegoInnerR(Float_t rad = 0.5)
	#tdrStyle.SetNumberContours(Int_t number = 20)

	tdrStyle.SetEndErrorSize(2)
	#  tdrStyle.SetErrorMarker(20)
	#tdrStyle.SetErrorX(0.)

	tdrStyle.SetMarkerStyle(20)

	#For the fit/function:
	tdrStyle.SetOptFit(1)
	tdrStyle.SetFitFormat("5.4g")
	tdrStyle.SetFuncColor(2)
	tdrStyle.SetFuncStyle(1)
	tdrStyle.SetFuncWidth(1)

	#For the date:
	tdrStyle.SetOptDate(0)
	#tdrStyle.SetDateX(Float_t x = 0.01)
	#tdrStyle.SetDateY(Float_t y = 0.01)

	# For the statistics box:
	tdrStyle.SetOptFile(0)
	tdrStyle.SetOptStat(0)  # To display the mean and RMS:   SetOptStat("mr")
	tdrStyle.SetStatColor(ROOT.kWhite)
	tdrStyle.SetStatFont(42)
	tdrStyle.SetStatFontSize(0.025)
	tdrStyle.SetStatTextColor(1)
	tdrStyle.SetStatFormat("6.4g")
	tdrStyle.SetStatBorderSize(1)
	tdrStyle.SetStatH(0.1)
	tdrStyle.SetStatW(0.15)
	#tdrStyle.SetStatStyle(Style_t style = 1001)
	#tdrStyle.SetStatX(Float_t x = 0)
	#tdrStyle.SetStatY(Float_t y = 0)

	# Margins:
	tdrStyle.SetPadTopMargin(0.05)
	tdrStyle.SetPadBottomMargin(0.13)
	tdrStyle.SetPadLeftMargin(0.16)
	tdrStyle.SetPadRightMargin(0.03)

	# For the Global title:
	tdrStyle.SetOptTitle(0)
	tdrStyle.SetTitleFont(42)
	tdrStyle.SetTitleColor(1)
	tdrStyle.SetTitleTextColor(1)
	tdrStyle.SetTitleFillColor(10)
	tdrStyle.SetTitleFontSize(0.05)
	#tdrStyle.SetTitleH(0)  # Set the height of the title box
	#tdrStyle.SetTitleW(0)  # Set the width of the title box
	#tdrStyle.SetTitleX(0)  # Set the position of the title box
	#tdrStyle.SetTitleY(0.985)  # Set the position of the title box
	#tdrStyle.SetTitleStyle(Style_t style = 1001)
	#tdrStyle.SetTitleBorderSize(2)

	# For the axis titles:
	tdrStyle.SetTitleColor(1, "XYZ")
	tdrStyle.SetTitleFont(42, "XYZ")
	tdrStyle.SetTitleSize(0.06, "XYZ")
	#tdrStyle.SetTitleXSize(Float_t size = 0.02)  # Another way to set the size?
	#tdrStyle.SetTitleYSize(Float_t size = 0.02)
	tdrStyle.SetTitleXOffset(0.9)
	tdrStyle.SetTitleYOffset(1.25)
	#tdrStyle.SetTitleOffset(1.1, "Y")  # Another way to set the Offset

	# For the axis labels:
	tdrStyle.SetLabelColor(1, "XYZ")
	tdrStyle.SetLabelFont(42, "XYZ")
	tdrStyle.SetLabelOffset(0.007, "XYZ")
	tdrStyle.SetLabelSize(0.05, "XYZ")

	# For the axis:
	tdrStyle.SetAxisColor(1, "XYZ")
	tdrStyle.SetStripDecimals(ROOT.kTRUE)
	tdrStyle.SetTickLength(0.03, "XYZ")
	tdrStyle.SetNdivisions(510, "XYZ")
	tdrStyle.SetPadTickX(1)   # To get tick marks on the opposite side of the frame
	tdrStyle.SetPadTickY(1)

	# Change for log plots:
	tdrStyle.SetOptLogx(0)
	tdrStyle.SetOptLogy(1)
	tdrStyle.SetOptLogz(0)

	# Postscript options:
	tdrStyle.SetPaperSize(20.,20.)
	#tdrStyle.SetLineScalePS(Float_t scale = 3)
	#tdrStyle.SetLineStyleString(Int_t i, const char* text)
	#tdrStyle.SetHeaderPS(const char* header)
	#tdrStyle.SetTitlePS(const char* pstitle)

	#tdrStyle.SetBarOffset(Float_t baroff = 0.5)
	#tdrStyle.SetBarWidth(Float_t barwidth = 0.5)
	#tdrStyle.SetPaintTextFormat(const char* format = "g")
	tdrStyle.SetPalette(56)
	tdrStyle.SetNumberContours(999)
	#tdrStyle.SetTimeOffset(Double_t toffset)
	#tdrStyle.SetHistMinimumZero(kTRUE)

	tdrStyle.cd()

	ROOT.gROOT.SetBatch()
	return tdrStyle


def tdrStyle2D():
	style = tdrStyle()
	style.SetOptLogy(0)
	style.SetOptLogz(1)
	#style.SetPadBottomMargin(0.09)
	style.SetPadTopMargin(0.06)
	style.SetPadRightMargin(0.18)
	#style.SetPadLeftMargin(0.09)
	style.SetCanvasDefW(1000)  #Width of canvas

	style.cd()
	return style

def setPalette(name="palette", ncontours=999):

	from array import array
	from ROOT import TColor, gStyle

	"""Set a color palette from a given RGB list
	stops, red, green and blue should all be lists of the same length
	see set_decent_colors for an example"""

	if name == "gray" or name == "grayscale":
		stops = [0.00, 0.34, 0.61, 0.84, 1.00]
		red   = [1.00, 0.84, 0.61, 0.34, 0.00]
		green = [1.00, 0.84, 0.61, 0.34, 0.00]
		blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
	elif name == "red":
		stops = [0, 1]
		red   = [1, 1]
		green = [1, 0]
		blue  = [0, 0]
	else:
		# default palette, looks cool
		stops = [0.00, 0.34, 0.61, 0.84, 1.00]
		red   = [0.00, 0.00, 0.87, 1.00, 0.51]
		green = [0.00, 0.81, 1.00, 0.20, 0.00]
		blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

	s = array('d', stops)
	r = array('d', red)
	g = array('d', green)
	b = array('d', blue)

	npoints = len(s)
	TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)

	gStyle.SetNumberContours(ncontours)


