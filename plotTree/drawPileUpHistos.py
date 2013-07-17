from treeFunctions import readHisto
import multiplot
import ROOT
import Styles
Styles.tdrStyle()
ROOT.gROOT.SetBatch()

data = readHisto( "../TreeWriter/pileUpReweighting/PU_dist.root" ,"pileup" )
data.SetLineColor(1)
data.SetMarkerColor(1)
data.SetMarkerStyle(20)
data.SetMarkerSize(.8)
data.Sumw2()

mc = readHisto( "../TreeWriter/pileUpReweighting/mc_PU_dist_S10.root" ,"pileup" )
mc.SetLineColor(2)
mc.SetMarkerColor(2)

for h in [mc, data ]:
	h.Scale( 1./h.Integral() )
	h.SetLineWidth(2)
	h.SetTitle( ";number of vertices;normalized entries" )
	h.SetTitleOffset( 1.5, "Y" )


muhist = multiplot.Multihisto()
muhist.addHisto( data, "Data", draw="ep" )
muhist.addHisto( mc, "MC" )

can = ROOT.TCanvas()
can.cd()
can.SetLogy(0)
muhist.Draw()
can.SaveAs("plots/puDistribution.pdf")
