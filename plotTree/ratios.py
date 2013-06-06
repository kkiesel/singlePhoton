#=======================================================
# Project: LocalAnalysis
#              SUSY Same Sign Dilepton Analysis
#
# File: ratios.py
#
# Author: Daniel Sprenger
#         daniel.sprenger@cern.ch
#=======================================================

import ROOT
import numpy as np

from array import array
import math



class Ratio:
    def __init__(self, numerator, denominator, numeratorSquaredError, denominatorSquaredError, xPos, width):
        self.numerators = [numerator]
        self.denominators = [denominator]
        self.numeratorSquaredErrors = [numeratorSquaredError]
        self.denominatorSquaredErrors = [denominatorSquaredError]
        self.xPos = [xPos]
        self.width = width

    @property
    def isValid(self):
        return (self.sumDenominator > 0)

    @property
    def ratio(self):
        value = -1.0
        if (self.isValid):
            value = self.sumNumerator / self.sumDenominator
        return value

    @property
    def sumNumerator(self):
        value = 0.0
        for num in self.numerators:
            value += num
        return value

    @property
    def sumDenominator(self):
        value = 0.0
        for den in self.denominators:
            value += den
        return value

    @property
    def sumNumeratorSquaredErrors(self):
        value = 0.0
        for numError in self.numeratorSquaredErrors:
            value += numError
        return value

    @property
    def sumDenominatorSquaredErrors(self):
        value = 0.0
        for denError in self.denominatorSquaredErrors:
            value += denError
        return value

    @property
    def errorX(self):
        return 0.5 * self.width

    @property
    def errorY(self):
        value = 0.0
        if (self.isValid and self.sumNumerator > 0):
            value = self.ratio * math.sqrt(self.sumNumeratorSquaredErrors / math.pow(self.sumNumerator, 2.0) + self.sumDenominatorSquaredErrors / math.pow(self.sumDenominator, 2.0))
        return value

    @property
    def xCenter(self):
        value = 0.0
        for x in self.xPos:
            value += x / len(self.xPos)
        return value

    def isFullEnough(self, rebinErrorBoundary=0.20):
        #if (self.sumNumerator > 0):
        return (self.sumNumerator > 0 and math.sqrt(self.sumNumeratorSquaredErrors) / self.sumNumerator <= rebinErrorBoundary)

    def addRatio(self, ratio):
        self.numerators.extend(ratio.numerators)
        self.denominators.extend(ratio.denominators)
        self.numeratorSquaredErrors.extend(ratio.numeratorSquaredErrors)
        self.denominatorSquaredErrors.extend(ratio.denominatorSquaredErrors)
        self.xPos.extend(ratio.xPos)
        self.width += ratio.width


class RatioGraph:
    def __init__(self, numerator, denominator, xMin=None, xMax=None):
        self.denominator = denominator
        self.numerator = numerator
        if xMin == None:
            self.xMin = numerator.GetBinLowEdge(1)
        else:
            self.xMin = xMin
        if xMax == None:
            self.xMax = numerator.GetBinLowEdge( numerator.GetNbinsX()+1 )
        else:
            self.xMax = xMax
        self.errors = []
        return

    def addErrorBySize(self, name, size, color=None, fillStyle=None, add=True):
        error = RatioError(name, self.xMin, self.xMax, size=size, add=add)
        if (color != None):
            error.color = color
        if (fillStyle != None):
            error.fillStyle = fillStyle
        self.errors.append(error)

    def addErrorByHistograms(self, name, denominatorUp, denominatorDown, color=None, fillStyle=None):
        error = RatioError(name, self.xMin, self.xMax, denominator=self.denominator, denominatorUp=denominatorUp, denominatorDown=denominatorDown)
        if (color != None):
            error.color = color
        if (fillStyle != None):
            error.fillStyle = fillStyle
        self.errors.append(error)

    def getGraph(self, adaptiveBinning=True, errors="yx"):
        ratios = []

        tempRatio = None
        for iBin in range(1, 1 + self.numerator.GetNbinsX()):
            num = self.numerator.GetBinContent(iBin)
            numError = self.numerator.GetBinError(iBin)
            den = self.denominator.GetBinContent(iBin)
            denError = self.denominator.GetBinError(iBin)
            x = self.numerator.GetBinCenter(iBin)
            width = self.numerator.GetBinWidth(iBin)

            # assure that bin is in view range
            # ignore empty starting bins
            if (self.xMin < x and x < self.xMax
                    and not (len(ratios) == 0 and tempRatio == None and den == 0)):

                if (tempRatio != None):
                    ratio = Ratio(num, den, math.pow(numError, 2.0), math.pow(denError, 2.0), x, width)
                    tempRatio.addRatio(ratio)
                else:
                    tempRatio = Ratio(num, den, math.pow(numError, 2.0), math.pow(denError, 2.0), x, width)

                # Start with the next bin if adaptive binning is disabled or the bin is full enough
                if not adaptiveBinning or (tempRatio.isFullEnough()):
                    ratios.append(tempRatio)
                    tempRatio = None

        if (tempRatio != None):
            ratios.append(tempRatio)

        xs = []
        ys = []
        yErrors = []
        widths = []
        for ratio in ratios:
            xs.append(ratio.xCenter)
            ys.append(ratio.ratio)
            yErrors.append(ratio.errorY)
            widths.append(ratio.errorX)



        # what did I mean with this?
        if "x" in errors and "y" in errors:
            graph = ROOT.TGraphAsymmErrors(len(xs), array("d", xs), array("d", ys), array("d", widths), array("d", widths), array("d", yErrors), array("d", yErrors))
        elif "x" in errors:
            graph = ROOT.TGraphAsymmErrors(len(xs), array("d", xs), array("d", ys), array("d", widths), array("d", widths), np.zeros(len(xs)), np.zeros(len(xs)))
        elif "y" in errors:
            graph = ROOT.TGraphAsymmErrors(len(xs), array("d", xs), array("d", ys), np.zeros(len(xs)), np.zeros(len(xs)), array("d", yErrors), array("d", yErrors))
        else:
            graph = ROOT.TGraphAsymmErrors(len(xs), array("d", xs), array("d", ys), np.zeros(len(xs)), np.zeros(len(xs)), array("d", 0), array("d", 0))

        return graph

    def getErrorGraphs(self):
        errorGraphs = []
        for iError, error in enumerate(self.errors):
            if (error.size != None):
                if (error.add):
                    print("Ignoring error, assuming it has already been added")
                else:
                    print("Adding error '%s' with size %f" % (error.name, error.size))
                    xCenter = 0.5 * (error.xMin + error.xMax)
                    xWidth = 0.5 * (error.xMax - error.xMin)
                    graph = ROOT.TGraphErrors(1, array("d", [xCenter]), array("d", [1.0]), array("d", [xWidth]), array("d", [error.size]))
                    graph.SetFillColor(error.color)
                    errorGraphs.append(graph)
            elif (error.hasHistograms):
                errorsUp = error.errorsUp
                errorsDown = error.errorsDown

                xs = []
                ys = []
                widths = []
                upErrors = []
                downErrors = []
                for (errorUp, errorDown) in zip(errorsUp, errorsDown):
                    if (errorUp.ratio >= 0.0 and errorDown.ratio >= 0.0):
                        xs.append(errorUp.xCenter)
                        ys.append(1.0)
                        widths.append(errorUp.errorX)

                        if (errorUp.ratio > 1.0):
                            if (errorDown.ratio > 1.0):
                                upErrors.append(max(errorUp.ratio - 1.0, errorDown.ratio - 1.0))
                                downErrors.append(0.0)
                            else:
                                upErrors.append(errorUp.ratio - 1.0)
                                downErrors.append(1.0 - errorDown.ratio)
                        else:
                            if (errorDown.ratio > 1.0):
                                upErrors.append(errorDown.ratio - 1.0)
                                downErrors.append(1.0 - errorUp.ratio)
                            else:
                                upErrors.append(0.0)
                                downErrors.append(max(1.0 - errorUp.ratio, 1.0 - errorDown.ratio))

                if (iError + 1 < len(self.errors) and self.errors[iError + 1].add):
                    if (self.errors[iError + 1].size != None):
                        size = self.errors[iError + 1].size
                        upErrors = [math.sqrt(prev ** 2 + size ** 2) for prev in upErrors]
                        downErrors = [math.sqrt(prev ** 2 + size ** 2) for prev in downErrors]
                    else:
                        print("Uncertainty to be added does not have fixed size. Adding not implemented, yet.")

                graph = ROOT.TGraphAsymmErrors(len(xs), array("d", xs), array("d", ys), array("d", widths), array("d", widths), array("d", downErrors), array("d", upErrors))
                graph.SetFillColor(error.color)
                graph.SetFillStyle(error.fillStyle)
                errorGraphs.append(graph)

        return errorGraphs

    def draw(self, pad, yMin=0.0, yMax=2.0, adaptiveBinning=True, errors="yx"):
        pad.cd()

        # axis
        nBinsX = 20
        nBinsY = 10
        self.hAxis = ROOT.TH2F("hAxis", "", nBinsX, self.xMin, self.xMax, nBinsY, yMin, yMax)
        self.hAxis.Draw("AXIS")
        self.hAxis.GetYaxis().SetNdivisions(int(yMax), 0, 2)
        self.hAxis.SetTitleOffset(0.3, "Y")
        self.hAxis.SetTitleSize(0.15, "Y")
        self.hAxis.SetYTitle("Ratio")
        self.hAxis.GetXaxis().SetLabelSize(0.0)
        self.hAxis.GetYaxis().SetLabelSize(0.15)

        self.errorGraphs = self.getErrorGraphs()
        self.errorGraphs.reverse()
        for errorGraph in self.errorGraphs:
            errorGraph.Draw("SAME2")

        self.oneLine = ROOT.TLine(self.xMin, 1.0, self.xMax, 1.0)
        self.oneLine.SetLineStyle(2)
        self.oneLine.Draw()

        self.hAxis.Draw("SAMEAXIS")

        self.graph = self.getGraph(adaptiveBinning, errors)
        self.graph.Draw("SAMEZ")

        pad.Update()


class RatioError:
    def __init__(self, name, xMin, xMax, size=None, denominator=None, denominatorUp=None, denominatorDown=None, add=False):
        self.name = name

        self.denominator = denominator
        self.denominatorUp = denominatorUp
        self.denominatorDown = denominatorDown
        self.size = size
        self.xMin = xMin
        self.xMax = xMax
        self.rebinErrorBoundary = 0.05

        self.add = add

        self.__ratiosUp__ = []
        self.__ratiosDown__ = []
        self.__errorsUp__ = None
        self.__errorsDown__ = None

        self.color = myColors['DarkSlateBlue']
        self.fillStyle = 1001

    @property
    def errorsUp(self):
        if (self.__errorsUp__ == None):
            if (len(self.__ratiosUp__) == 0):
                self._calculateRatios()

            self.__errorsUp__ = []
            for ratio in self.__ratiosUp__:
                self.__errorsUp__.append(ratio)

        return self.__errorsUp__

    @property
    def errorsDown(self):
        if (self.__errorsDown__ == None):
            if (len(self.__ratiosDown__) == 0):
                self._calculateRatios()

            self.__errorsDown__ = []
            for ratio in self.__ratiosDown__:
                self.__errorsDown__.append(ratio)

        return self.__errorsDown__

    @property
    def hasHistograms(self):
        return (self.denominator != None and self.denominatorUp != None and self.denominatorDown != None)

    def _calculateRatios(self):
        if (self.hasHistograms):
            tempRatioUp = None
            tempRatioDown = None
            for iBin in range(1, 1 + self.denominator.GetNbinsX()):
                den = self.denominator.GetBinContent(iBin)
                denError = self.denominator.GetBinError(iBin)
                denUp = self.denominatorUp.GetBinContent(iBin)
                denUpError = self.denominatorUp.GetBinError(iBin)
                denDown = self.denominatorDown.GetBinContent(iBin)
                denDownError = self.denominatorDown.GetBinError(iBin)


                x = self.denominator.GetBinCenter(iBin)
                width = self.denominator.GetBinWidth(iBin)

                # assure that bin is in view range
                if (self.xMin < x and x < self.xMax):
                    ratioUp = Ratio(denUp, den, math.pow(denUpError, 2.0), math.pow(denError, 2.0), x, width)
                    ratioDown = Ratio(denDown, den, math.pow(denDownError, 2.0), math.pow(denError, 2.0), x, width)


                    if (tempRatioUp != None):
                        tempRatioUp.addRatio(ratioUp)
                        tempRatioDown.addRatio(ratioDown)
                    else:
                        tempRatioUp = ratioUp
                        tempRatioDown = ratioDown

                    if (tempRatioUp.isFullEnough(self.rebinErrorBoundary) and tempRatioDown.isFullEnough(self.rebinErrorBoundary)):
                        self.__ratiosUp__.append(tempRatioUp)
                        self.__ratiosDown__.append(tempRatioDown)
                        tempRatioUp = None
                        tempRatioDown = None

            if (tempRatioUp != None):
                self.__ratiosUp__.append(tempRatioUp)
                self.__ratiosDown__.append(tempRatioDown)
        else:
            print("Trying to calculate error ratios, but histograms not set!")



