#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from treeFunctions import *

filename = "W_1700_720_375_V01.cutFlow_tree.root"

can = ROOT.TCanvas( randomName(), "", 600, 600 )
can.SetLogy(0)
can.SetRightMargin( 0.08 )
can.SetLeftMargin( 0.12 )
can.SetBottomMargin( 0.08 )

h = readHisto( filename, "cutFlow" )

# scale
#lumi = 19712
#xsec = 0.3164
nGen = 60000
#h.Scale( lumi*xsec/nGen )

#triggerEfficiency = 0.972
#triggerEfficiency = 1
#h.Scale( 1./triggerEfficiency )

#prettify labels
labels = {
	"lumiSection": "",
	"vertex": "n_{vertex} #geq 1",
	"photon": "n_{#gamma} #geq 1",
	"no e": "n_{e} = 0",
	"no #mu": "n_{#mu} = 0",
	"nJets": "n_{jets} #geq 2",
	"H_{T}": "H_{T}#geq500 GeV",
	"#slash{E}_{T}": "#slash{E}_{T}#geq100 GeV"
}
for bin in range(h.GetNbinsX()):
	oldLabel = h.GetXaxis().GetBinLabel( bin )
	if oldLabel in labels:
		h.GetXaxis().SetBinLabel( bin, labels[oldLabel] )

# only interesting range
h.GetXaxis().SetRange(2,9)

# titel and ylabel
h.GetYaxis().SetTitle( "Events" )
h.GetYaxis().SetTitleOffset(1.6)

h.SetLineColor(1)
h.SetLineWidth(2)
h.Draw("hist")

infoText = ROOT.TLatex()
infoText.SetNDC()
infoText.SetTextFont( h.GetLabelFont() )
infoText.SetTextSize( h.GetLabelSize()*0.8 )
infoText.SetText( .01, .97, "CMS Simulation: GGM Wino-like #tilde{#chi}_{1}^{0}          19.7fb^{-1} (8 TeV)" )
infoText.Draw()

infoText2 = ROOT.TLatex()
infoText2.SetNDC()
infoText2.SetTextFont( h.GetLabelFont() )
infoText2.SetTextSize( h.GetLabelSize()*0.7 )
infoText2.SetText( .7, .8, "particle  mass [GeV]" )
infoText2.Draw()
lineSpacing=0.0
for text1, text2 in [ ("#tilde{q} ", "1700"),("#tilde{g} "," 720"), ("#tilde{#chi}^{0}_{1}", " 375")]:
	lineSpacing -= 0.04
	infoText2.DrawLatex( infoText2.GetX(), infoText2.GetY()+lineSpacing, text1+"           "+text2 )

ROOT.gPad.SaveAs("plots/cutFlow%s.pdf"%filename.replace(".root", "") )

for bin in range(h.GetNbinsX()):
	print h.GetXaxis().GetBinLabel(bin), "%.2f"%(h.GetBinContent(bin)/nGen*100), "%.2f"%(sqrt(h.GetBinContent(bin))/nGen*100)
