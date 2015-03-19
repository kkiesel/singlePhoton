import ROOT

class Multiplot:
    def __init__( self ):
        self.hists = []
        self.histsToStack = []

        # todo: impmlement setter and getter
        self.minimum = None
        self.maximum = None

        self.leg = ROOT.TLegend(.7,.7,.93,.92)

    def add( self, h ):
        # important histograms first, they will be drawn on top
        self.hists.append( h )

    def addStack( self, h ):
        self.histsToStack.append( h )

    def getMinimum( self ):
        return min( [ h.GetMinimum(0) for h in self.hists ] )

    def getMaximum( self ):
        return max( [ h.GetMaximum(0) for h in self.hists ] )

    def stackHists( self ):
        if not self.histsToStack:
            return
        #stack = ROOT.THStack()
        #stack.SetTitle( ";%s;%s"%(self.histosToStack[0][0].GetXaxis().GetTitle(),self.histosToStack[0][0].GetYaxis().GetTitle()) )
        stack = ROOT.THStack( self.histsToStack[0] ) # get title, etc
        for h in histsToStack:
            h.SetFillColor( h.GetLineColor() )
            h.SetLineColor( ROOT.kBlack )
            stack.Add( h )

        self.hists.append( stack )


    def Draw( self ):
        self.stackHists()

        minimum = self.getMinimum()
        maximum = self.getMaximum()
        if not ROOT.gPad or ROOT.gPad.GetLogy():
            maximum = 2.5*maximum
            minimum = 0.3*minimum
        else:
            maximum = maximum + (maximum-minimum)*.1
            minimum = minimum - (maximum-minimum)*.1

        if self.maximum != None:
            maximum = self.maximum
        if self.minimum != None:
            minimum = self.minimum

        # fill legend (in correct order)
        for h in self.hists:
            # todo: special cases?
            self.leg.AddEntry( h, h.GetTitle(), "lpe" )
        for h in self.histsToStack:
            self.leg.AddEntry( h, h.GetTitle(), "f" )

        # change the order for drawing
        h.reverse()
        hists[0].SetMinimum( minimum )
        hists[0].SetMaximum( maximum )
        hists[0].Draw()

        for h in hists[1:]:
            h.Draw("same")
