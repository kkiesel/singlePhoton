#!/usr/bin/env python2

import ROOT
import style
import math

metStr = "#it{E}_{#kern[-0.1]{T}} #scale[0.7]{#kern[-.5]{#lower[-0.5]{miss}}}"
#metStr = "#it{E}_{#kern[-0.1]{T}}^{miss}"

def bottomPad( relativeSize=.2 ):

    # Delete label and title of all histograms in the current pad
    for primitive in ROOT.gPad.GetListOfPrimitives():
        if getattr( primitive, "GetXaxis", None ):
            xaxis = primitive.GetXaxis()
            xaxis.SetLabelSize(0)
            xaxis.SetLabelColor(0)
            xaxis.SetLabelOffset(1000)
            xaxis.SetTitle("")
            xaxis.SetTitleColor(0)
            xaxis.SetTitleSize(0)

    ROOT.gPad.SetBottomMargin( relativeSize + (1-relativeSize)*ROOT.gPad.GetBottomMargin() - relativeSize*ROOT.gPad.GetTopMargin() )
    rPad = ROOT.TPad( "rPad", "ratio", 0, 0, 1, 1 )
    rPad.SetTopMargin( (1-relativeSize) - (1-relativeSize)*rPad.GetBottomMargin() + relativeSize*rPad.GetTopMargin() )
    rPad.SetFillStyle(3955)
    rPad.Draw()
    rPad.cd()
    rPad.SetLogy(0)

    return rPad

def horizontalLine( xaxis, pos=1 ):
    oneLine = ROOT.TLine()
    oneLine.SetLineStyle(2)
    xMin = xaxis.GetBinLowEdge(1)
    xMax = xaxis.GetBinUpEdge( xaxis.GetNbins() )
    oneLine.DrawLine( xMin, pos, xMax, pos )

def fillErrorHist( hist_up, hist_dn, newname="" ):
    if not newname: newname = hist_up.GetName()+"new"

    h = hist_up.Clone( newname )
    for bin in range( h.GetNbinsX()+2 ):
        up = hist_up.GetBinContent(bin)
        dn = hist_dn.GetBinContent(bin)
        h.SetBinContent( bin, 0.5*(up+dn) )
        h.SetBinError( bin, 0.5*(up-dn) )
    return h

def addUncertQuadratic( h1, h2, newname="" ):
    if not newname: newname = h1.GetName()+"new"

    h = h1.Clone( newname )
    for bin in range( h.GetNbinsX()+2 ):
        c1 = h1.GetBinContent(bin)
        c2 = h2.GetBinContent(bin)
        e1 = h1.GetBinError(bin)
        e2 = h2.GetBinError(bin)
        if c1 != c1: print "Warning: the content should be the same", c1, c2
        h.SetBinError( bin, math.sqrt(e1**2+e2**2) )
    return h

def subUncertQuadratic( h1, h2, newname="" ):
    if not newname: newname = h1.GetName()+"new"

    h = h1.Clone( newname )
    for bin in range( h.GetNbinsX()+2 ):
        c1 = h1.GetBinContent(bin)
        c2 = h2.GetBinContent(bin)
        e1 = h1.GetBinError(bin)
        e2 = h2.GetBinError(bin)
        if c1 != c1: print "Warning: the content should be the same", c1, c2
        h.SetBinError( bin, math.sqrt(e1**2-e2**2) )
    return h



###############################################################################
# Architecture:
#                       BeautyPlot
#         RazorPlot                       SingleGPlot
#   Signal    Closure                 Final      ewkclosure qcdClosure
# T5gg Bino  1 2 3 4 5                FinalPlot


class BeautyPlot:
    isSimulation = True
    selectionAppendix = ""

    labelTemplate = "#font[61]{CMS}"
    simulationStr = " #scale[0.76]{#font[52]{Simulation}}"
    lumiInfoText = "19.7 fb^{-1} (8 TeV)"
    cutLabex = ""

    def drawLabel( self ):
        label = ROOT.TLatex()
        if self.isSimulation: self.labelTemplate += self.simulationStr
        #self.selection += self.selectionAppendix
        label.DrawLatexNDC( ROOT.gPad.GetLeftMargin()+0.02, .895, self.labelTemplate )
        label.DrawLatexNDC( .68, .955, self.lumiInfoText )
        label.DrawLatexNDC( .68, .88, self.selection )
        if self.selectionAppendix:
            label.DrawLatexNDC( .68, .84, self.selectionAppendix )

    def readC( self ):
        ROOT.gROOT.ProcessLine( " .x %s"% self.infile )

    def save( self ):
        for ending in [ "pdf", "C" ]:
            ROOT.gPad.GetCanvas().SaveAs( "output/{}.{}".format( self.output, ending ) )

    def __init__( self ):
        self.readC()
        if ROOT.gDirectory:
            for obj in ROOT.gDirectory.GetList():
                if hasattr( obj, "UseCurrentStyle" ):
                    obj.UseCurrentStyle()

        c = ROOT.TCanvas()
        self.draw()
        self.drawLabel()
        self.save()
        if ROOT.gDirectory:
            ROOT.gDirectory.Clear()

############### Razor Ploting #################

class RazorPlot( BeautyPlot ):
    selection = "Razor #gamma#gamma #geq1jet"
    bottomPlotRange = -2.2, 2.2

class SignalComparison( RazorPlot ):


    def draw( self ):
        bkg = ROOT.hist_high
        bkg.SetMaximum(1e3)
        bkg.SetYTitle("Events")
        bkg.UseCurrentStyle()
        bkg.SetLineWidth(3)
        bkg.SetLineColor( ROOT.kBlack )
        s1 = ROOT.hist_high_sig1
        s1.SetLineColor( ROOT.kRed )
        s2 = ROOT.hist_high_sig2
        s2.SetLineColor( ROOT.kGreen )
        s3 = ROOT.hist_high_sig3
        s3.SetLineColor( ROOT.kBlue )

        ROOT.gROOT.cd()
        self.leg = ROOT.TLegend(.25,.59,.95,.87)
        self.leg.SetFillColor(0)
        self.leg.AddEntry( bkg, "Background model" )
        self.leg.AddEntry( 0, self.legEntries[0], "" )
        bkg.Draw()
        for i, h in enumerate( [s1, s2, s3] ):
            h.SetLineWidth(3)
            self.leg.AddEntry( h, self.legEntries[i+1], "l" )
            h.Draw("same")
        self.leg.Draw()

class T5SignalPlot( SignalComparison ):
    infile = "input/t5gg_s_v_b.C"
    output = "DiPhoton_SignalShape_T5gg"
    legEntries = [
        "T5gg signals:",
        "m_{#tilde{g}} = 1350 GeV, m_{#tilde{#chi}^{0}_{1}} =   225 GeV",
        "m_{#tilde{g}} = 1350 GeV, m_{#tilde{#chi}^{0}_{1}} =   675 GeV",
        "m_{#tilde{g}} = 1350 GeV, m_{#tilde{#chi}^{0}_{1}} = 1275 GeV"
    ]

class BinoSignalPlot( SignalComparison ):
    infile = "input/ra3_s_v_b.C"
    output = "DiPhoton_SignalShape_GGMbino"
    legEntries = [
        "GGMbino signals:",
        "m_{#tilde{q}} = 1500 GeV, m_{#tilde{g}} = 1820 GeV",
        "m_{#tilde{q}} = 1700 GeV, m_{#tilde{g}} = 1520 GeV",
        "m_{#tilde{q}} = 1900 GeV, m_{#tilde{g}} = 1320 GeV"
    ]

class RazorComparison( RazorPlot ):
    isSimulation = False
    def draw( self ):
        self.getPlots()
        self.pred.UseCurrentStyle()

        if False:
            # check z-score
            bin = 8
            print "Data", self.data.GetBinContent(bin+1)
            x,y, ey = ROOT.Double(), ROOT.Double(), ROOT.Double()
            self.predUncert.GetPoint( bin, x, y )
            print "pred", y, "m", self.predUncert.GetErrorYlow(bin), "p", self.predUncert.GetErrorYhigh(bin)
            self.ratio.GetPoint( bin, x, y )
            print "z-score", y
            print "z-score uncert", self.ratioSys.GetErrorYlow(bin)
            print "z-score uncert", self.ratioSys.GetErrorYhigh(bin)


        for h in self.data, self.ratio:
            h.SetMarkerStyle(20)
            h.SetMarkerSize(1)
            h.SetMarkerColor( ROOT.kBlack )
            h.SetLineColor( ROOT.kBlack )

        self.predUncert.SetMaximum(1e4)
        self.predUncert.SetMarkerStyle(24)
        self.predUncert.SetMarkerColor( ROOT.kBlack )

        self.data.SetBinErrorOption( 1 )
        self.pred.SetYTitle("Events")
        self.pred.Draw("p")
        self.predUncert.Draw("e2 p")

        if isinstance( self, RazorSignalInj ):
            self.signal.SetLineWidth(3)
            self.signal.SetLineColor( ROOT.kRed )
            self.signal.Draw("same")

        #self.pred.Draw("p same")
        self.data.Draw("e0 x0 p same")

        self.leg = ROOT.TLegend(.35,.65,.9,.85)
        self.leg.SetFillColor(0)
        self.leg.AddEntry( self.data, self.legEntries[0], "ep" )
        self.leg.AddEntry( self.predUncert, self.legEntries[1], "fp" )
        if isinstance( self, RazorSignalInj ):
            self.leg.SetY1(0.57)
            self.leg.AddEntry( self.signal, "GGMbino #scale[0.8]{#splitline{m_{#tilde{q}}=1400 GeV}{m_{#tilde{g}}=1820 GeV}}", "l" )
        self.leg.SetTextSize( ROOT.gStyle.GetTitleSize())
        self.leg.Draw()

        self.ratio.SetHistogram( self.data.Clone("ratioHist") )
        bottomPad()

        self.ratio.GetYaxis().SetNdivisions( 3, 0, 0, True )
        self.ratio.SetMaximum(2.5)
        self.ratio.SetMinimum(-2.5)
        if isinstance( self, RazorSignalInj ):
            self.ratio.SetMaximum(10)
            self.ratio.SetMinimum(-2)
            self.ratio.GetYaxis().SetNdivisions( 5, 0, 0, True )

        self.ratio.SetTitle(";M_{R} (TeV);z-score")
        self.ratio.SetLineColor( 0 )
        self.ratio.Draw()
        self.ratioSys.Draw("same e2")
        horizontalLine( self.ratio.GetXaxis(), 0)
        self.ratio.Draw("p same")


class RazorControlLowR( RazorComparison ):
    infile = "input/fake_low_fit_edit.C"
    output = "DiPhoton_Control_LowR"

    legEntries = [
        "Low R^{2} Control Sample",
        "Fit Prediction"
    ]

    def getPlots( self ):
        self.data = ROOT.hist_low_data
        self.pred = ROOT.hist_low_pred
        self.predUncert = ROOT.shapeUncert
        self.ratio = ROOT.ratio
        self.ratioSys = ROOT.ratio_shape



class RazorControlHighR( RazorComparison ):
    infile = "input/fake_extrap_edit.C"
    output = "DiPhoton_Control_HighR"

    legEntries = [
        "High R^{2} Control Sample",
        "Estimation in High R^{2}"
    ]
    def getPlots( self ):
        self.data = ROOT.hist_high_data
        self.pred = ROOT.hist_high_pred
        self.predUncert = ROOT.shapeUncert
        self.ratio = ROOT.ratio
        self.ratioSys = ROOT.ratio_shape

class RazorLowR( RazorComparison ):
    infile = "input/low_fit_edit.C"
    output = "DiPhoton_LowR"

    legEntries = [
        "Data (Low R^{2})",
        "Fit Prediction"
    ]

    def getPlots( self ):
        self.data = ROOT.hist_low_data
        self.pred = ROOT.hist_low_pred
        self.predUncert = ROOT.shapeUncert
        self.ratio = ROOT.ratio
        self.ratioSys = ROOT.ratio_shape

class RazorHighR( RazorComparison ):
    infile = "input/extrap_edit.C"
    output = "DiPhoton_HighR"

    legEntries = [
        "Data (High R^{2})",
        "Estimation in High R^{2}"
    ]
    def getPlots( self ):
        self.data = ROOT.hist_high_data
        self.pred = ROOT.hist_high_pred
        self.predUncert = ROOT.shapeUncert
        self.ratio = ROOT.ratio
        self.ratioSys = ROOT.ratio_shape

class RazorSignalInj( RazorComparison ):
    infile = "input/signal_inj_edit.C"
    output = "DiPhoton_SignalInjection"

    bottomPlotRange = -2.2, 9

    legEntries = [
        "Control + Signal (High R^{2})",
        "Estimation in High R^{2}"
    ]
    def getPlots( self ):
        self.data = ROOT.hist_high_data
        self.pred = ROOT.hist_high_pred
        self.predUncert = ROOT.shapeUncert
        self.ratio = ROOT.ratio
        self.ratioSys = ROOT.ratio_shape
        self.signal = ROOT.hist_signal





############### Single photon

class SingleClosure( BeautyPlot ):
    selection = "#geq1#gamma #geq2jet"
    bottomPlotRange = 0, 2.2
    dataset = ""

    legendPos1 = .5,.63

    def drawLegend( self ):
        dx = 0.4
        dy = 0.25
        legendPos2 = self.legendPos1[0]+dx, self.legendPos1[1]+dy
        self.leg = ROOT.TLegend( *(self.legendPos1+legendPos2) )
        self.leg.SetFillStyle(0)
        if self.dataset:
            self.leg.SetHeader( self.dataset )
        self.leg.AddEntry( self.data, "Direct Simulation", "pe" )
        self.leg.AddEntry( self.total, "Prediction #pm #sigma_{total}", "lf" )
        self.leg.AddEntry( self.syst, "#pm #sigma_{syst}", "f" )
        self.leg.Draw()

    def drawStacked( self ):
        pass

    def draw( self ):
        self.getPlots()
        self.total.UseCurrentStyle()
        self.ratioTotal.UseCurrentStyle()

        if isinstance( self, NJetPlot ):
            for h in [ self.ratioTotal, self.total]:
                ax = h.GetXaxis()
                ax.SetNdivisions( 1, 0, 0, False )
                for bin in range(1, h.GetNbinsX()+1):
                    ax.SetBinLabel( bin, str(bin+1) )
                ax.SetLabelSize( h.GetYaxis().GetLabelSize()*1.5 )
                ax.SetLabelFont( h.GetYaxis().GetLabelFont() )
                ax.SetLabelOffset( 0.02 )
        if isinstance( self, MetPlot ):
            for h in [ self.ratioTotal, self.total]:
                ax = h.GetXaxis()
                ax.SetNdivisions( 8, 5, 0 )

        if isinstance( self, EwkClosureMet ):
            self.total.SetMinimum(1.1e-2)


        # set plot options
        for h in self.data, self.ratio:
            h.SetLineColor(1)
            h.SetLineWidth(3)
            h.SetMarkerStyle(20)

        for h in self.total, self.ratioTotal:
            h.SetFillStyle( 3245 )
            h.SetFillColor( ROOT.kBlue )
            h.SetLineWidth(0)
            h.SetLineColor( ROOT.kBlue )

        for h in self.syst, self.ratioSyst:
            h.SetFillStyle( 3254 )
            h.SetFillColor( ROOT.kRed )
            h.SetLineWidth(0)
            h.SetLineColor(0)

        for h in self.stat, self.ratioStat:
            h.SetMarkerSize(0)
            h.SetMarkerColor( ROOT.kRed )
            h.SetLineColor( ROOT.kRed )
            h.SetLineWidth(5)

        for h in self.pred, self.total:
            self.pred.SetLineColor( ROOT.kBlue )
            self.pred.SetLineWidth( 1 )

        self.pred.SetFillStyle( 0 )

        self.total.SetTitle( self.title )
        self.total.Draw("e2")
        self.drawStacked()
        self.total.Draw("e2 same")
        self.syst.Draw("e2 same")
        self.stat.Draw("e x0 same")
        self.pred.SetDrawOption("hist same")
        self.pred.Draw("hist same")
        self.data.Draw("e0p x0 same")

        self.drawLegend()

        self.leg2 = ROOT.TLegend( self.leg )
        self.leg2.SetX1( .5*(self.leg.GetX1()+self.leg.GetX2()) )
        self.leg2.Clear()
        self.leg2.AddEntry( 0, " ", "" )
        self.leg2.AddEntry( 0, " ", "" )
        if not isinstance( self, FinalPlotMet ):
            self.leg2.AddEntry( 0, " ", "" )
        self.leg2.AddEntry( self.stat, "#pm #sigma_{stat}", "e" )
        if isinstance( self, FinalPlotMet ):
            self.leg2.AddEntry( 0, " ", "" )
            self.leg2.AddEntry( 0, " ", "" )
            self.leg2.AddEntry( 0, " ", "" )
            self.leg2.AddEntry( 0, " ", "" )

        self.leg2.Draw()


        bottomPad()

        self.ratioTotalUnity = self.ratioTotal.Clone("ratioTotalUnity")
        self.ratioTotalUnity.SetFillStyle(0)

        self.ratioTotal.SetTitle( self.title )
        self.ratioTotal.GetYaxis().SetNdivisions( 3, 5, 0, True )
        self.ratioTotal.GetYaxis().SetRangeUser( 0, 2.1 )
        self.ratioTotal.GetYaxis().SetTitle( "Sim. / Pred." )
        if hasattr( self, "qcd" ):
            self.ratioTotal.GetYaxis().SetTitle( "Data / Pred." )
        self.ratioTotal.Draw("e2")
        self.ratioTotalUnity.Draw("hist same")
        self.ratioSyst.Draw("e2 same")
        self.ratioStat.Draw("e x0 same")
        #horizontalLine( self.ratio.GetXaxis(), 1)
        self.ratio.Draw("e0p same x0")

class MetPlot():
    title = ";%s (GeV);Events / GeV"%metStr

class NJetPlot():
    title = ";Jets;Events"
    legendPos1 = .2,.35

class EwkClosure():
    dataset = "t#bar{t}, W"

    def getPlots( self ):
        self.data = ROOT.gHist
        self.stat = ROOT.eHist.Clone("stat_eHist")
        self.syst = ROOT.eHistSys
        self.pred = ROOT.eHist
        self.total = addUncertQuadratic( self.syst, self.pred )

        # ratio
        self.ratio = ROOT.ratio_gHist
        self.ratioSyst = ROOT.ratioSys_eHist
        self.ratioTotal = ROOT.ratioTot_eHist
        self.ratioStat = subUncertQuadratic( self.ratioTotal, self.ratioSyst )

class QcdClosure():
    dataset = "Multijet, #gamma+jet"


class EwkClosureMet( EwkClosure, MetPlot, SingleClosure ):
    infile = "input/ewkClosure_TTJetsW_met.C"
    output = "SinglePhoton_EWKclosure_met"

class EwkClosureNjet( EwkClosure, NJetPlot, SingleClosure ):
    infile = "input/ewkClosure_TTJetsW_nGoodJets.C"
    output = "SinglePhoton_EWKclosure_nJet"



class QcdClosureMet( QcdClosure, MetPlot, SingleClosure ):
    infile = "input/Closure_Combined_met_log_edit.C"
    output = "SinglePhoton_QCDclosure_met"

    def getPlots( self ):

        # normal
        totalUp = ROOT.met_1090__3832
        systUp = ROOT.we_met_1091__3833
        systDn = ROOT.we_met_1091__3834
        totalDn = ROOT.met_1090__3835

        self.syst = fillErrorHist( systUp, systDn )
        self.stat = ROOT.hs.GetStack().Last()
        self.total = fillErrorHist( totalUp, totalDn )
        self.pred = self.total.Clone("pred")
        self.data = ROOT.met_1092__3836

        self.total.SetMaximum( self.total.GetMaximum()*10 )

        # ratio
        ratio_totalUp = ROOT.met_1090__3838
        ratio_systUp = ROOT.we_met_1091_1
        ratio_systDn = ROOT.we_met_1091
        ratio_totalDn = ROOT.met_1090__3839

        self.ratioSyst = fillErrorHist( ratio_systUp, ratio_systDn )
        self.ratioTotal = fillErrorHist( ratio_totalUp, ratio_totalDn )
        self.ratio = ROOT.met_1092
        self.ratioStat = ROOT.met_1090


class QcdClosureNjet( QcdClosure, NJetPlot, SingleClosure ):
    infile = "input/Final_Combined_n_jet7_log_edit.C"
    output = "SinglePhoton_QCDclosure_njet"
    selectionAppendix = "%s #geq 100GeV"%metStr

    def getPlots( self ):

        # normal
        totalUp = ROOT.n_jet7_1684__6802
        systUp = ROOT.we_n_jet7_1685__6803
        systDn = ROOT.we_n_jet7_1685__6804
        totalDn = ROOT.n_jet7_1684__6805

        self.syst = fillErrorHist( systUp, systDn )
        self.stat = ROOT.hs.GetStack().Last()
        self.total = fillErrorHist( totalUp, totalDn )
        self.pred = self.total.Clone("pred")
        self.data = ROOT.n_jet7_1686__6806

        self.total.SetMaximum( self.total.GetMaximum()*3 )

        # ratio
        ratio_totalUp = ROOT.n_jet7_1684__6808
        ratio_systUp = ROOT.we_n_jet7_1685_1
        ratio_systDn = ROOT.we_n_jet7_1685
        ratio_totalDn = ROOT.n_jet7_1684__6809

        self.ratioSyst = fillErrorHist( ratio_systUp, ratio_systDn )
        self.ratioTotal = fillErrorHist( ratio_totalUp, ratio_totalDn )
        self.ratio = ROOT.n_jet7_1686
        self.ratioStat = ROOT.n_jet7_1684

class FinalPlotMet( MetPlot, SingleClosure ):
    infile = "input/Closure_Data_met_log_edit.C"
    output = "SinglePhoton_Data_met"
    isSimulation = False

    def getPlots( self ):

        # normal
        totalUp = ROOT.met_495__1509
        systUp = ROOT.we_met_496__1510
        systDn = ROOT.we_met_496__1511
        totalDn = ROOT.met_495__1512

        self.syst = fillErrorHist( systUp, systDn )
        self.stat = ROOT.hs.GetStack().Last()
        self.total = fillErrorHist( totalUp, totalDn )
        self.pred = self.total.Clone("pred")
        self.data = ROOT.met_497__1514

        # extra plots
        self.signal = ROOT.met_498__1513
        self.isr = ROOT.met_499
        self.ewk = ROOT.met_501
        self.qcd = ROOT.met_495__1518
        self.qcd.Scale( 1, "width" )

        # ratio
        ratio_totalUp = ROOT.met_495__1516
        ratio_systUp = ROOT.we_met_496_1
        ratio_systDn = ROOT.we_met_496
        ratio_totalDn = ROOT.met_495__1517

        self.ratioSyst = fillErrorHist( ratio_systUp, ratio_systDn )
        self.ratioTotal = fillErrorHist( ratio_totalUp, ratio_totalDn )
        self.ratioStat = ROOT.met_495_2
        self.ratio = ROOT.met_497


    def drawStacked( self ):


        self.isr.SetLineColor( ROOT.kOrange + 7 )
        self.ewk.SetLineColor( ROOT.kGreen + 3 )
        self.qcd.SetLineColor( ROOT.kCyan - 10 )
        self.qcd.SetLineColor( ROOT.kAzure + 10 )

        for h in self.isr, self.ewk, self.qcd:
            h.SetFillStyle( 1001 )
            h.SetFillColor( h.GetLineColor() )
            h.SetLineWidth(1)
            h.SetLineColor(1)
            h.SetMarkerSize(0)
            h.Sumw2(0)

        self.signal.Sumw2(0)
        self.signal.SetLineColor( ROOT.kRed )
        self.signal.SetLineWidth(2)

        self.hs = ROOT.THStack()
        self.hs.UseCurrentStyle()
        self.hs.Add( self.ewk )
        self.hs.Add( self.isr )
        self.hs.Add( self.qcd )
        self.hs.Add( self.signal )
        self.hs.SetMinimum(0.06)
        self.hs.SetMaximum(1.5e4)
        self.hs.SetTitle( self.title )
        self.hs.Draw()

    def drawLegend( self ):
        self.leg = ROOT.TLegend(.47,.53,.95,.87)
        self.leg.SetFillStyle(0)
        if self.dataset:
            self.leg.SetHeader( self.dataset )
        self.leg.AddEntry( self.data, "Data", "pe" )
        self.leg.AddEntry( self.total, "Prediction #pm #sigma_{total}", "lf" )
        self.leg.AddEntry( self.syst, "#pm #sigma_{syst}", "f" )
        self.leg.AddEntry( self.qcd, "Multijet (+#gamma)", "f" )
        self.leg.AddEntry( self.isr, "Z#gamma, W#gamma, t#bar{t}#gamma", "f" )
        self.leg.AddEntry( self.ewk, "e#rightarrow#gamma", "f" )
        self.leg.AddEntry( self.signal, "GGMwino #scale[0.8]{#splitline{m_{#tilde{q}}=1700 GeV}{m_{#tilde{g}}= 720 GeV}}", "l" )
        self.leg.Draw()






if __name__ == "__main__":

    T5SignalPlot()
    BinoSignalPlot()

    RazorControlLowR()
    RazorControlHighR()

    RazorSignalInj()

    RazorLowR()
    RazorHighR()


    EwkClosureMet()
    EwkClosureNjet()
    QcdClosureMet()
    QcdClosureNjet()

    FinalPlotMet()

