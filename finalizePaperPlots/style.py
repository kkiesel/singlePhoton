import ROOT

s = ROOT.gStyle

s.SetCanvasDefH(600)
s.SetCanvasDefW(600)

s.SetOptLogy()

#s.SetPadTickX(1)
#s.SetPadTickY(1)
s.SetTickLength(0.02, "xyz")

s.SetPadTopMargin(0.06)
s.SetPadBottomMargin(0.10)
s.SetPadRightMargin(0.035)
s.SetPadLeftMargin(0.12)

s.SetHistTopMargin(0.9)

# text:
font = 42
size = 0.035*1.3

s.SetTextFont( font )
s.SetTitleFont( font, "xyz" )
s.SetLabelFont( font, "xyz" )
s.SetLegendFont( font )

s.SetTitleOffset( 1.3, "y" )

s.SetTextSize( size )
s.SetLabelSize( size, "xyz" )
s.SetTitleSize( size, "xyz" )


s.SetLegendBorderSize(0)
s.SetLegendFillColor(0)

ROOT.gROOT.SetBatch()
ROOT.gROOT.ForceStyle()


