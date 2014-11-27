#!/usr/bin/env python2
import ROOT
import re
import lhapdf
import pickle
import math

crossSectionFileBino = "../../infos/Spectra_gsq_B_8TeV.xsec"
crossSectionFileWino = "../../infos/Spectra_gsq_W_8TeV.xsec"
crossSectionFileT5 = "../../infos/simplifiedModel.xsec"
lumiData = 19712.

def readTree( filename, treename = "susyTree" ):
	"""
	filename: name of file containing the tree
	treename: name of the tree
	returns: TChain Object
	"""
	import ROOT
	tree = ROOT.TChain( treename )
	tree.AddFile( filename )
	return tree


def readSignalXSection( filename ):
    """Read xsection and other informations for various signal MC from a file
    found at https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA3PrivateSignalMC2012

    returns list[ (point1, point2) ] = ( sigmaNLO, errorUp, errorDown )
    """
    f = open( filename )
    text = f.readlines()
    f.close()

    info = {}

    floatMatch = "\-{0,1}\d\.\d+e[+-]\d{2}"
    import re
    for t in text:
        matches = re.match(" (\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+) LO: (%s) \+ (%s) \- (%s) NLO: (%s) \+ (%s) \- (%s)\n"%(floatMatch,floatMatch,floatMatch,floatMatch,floatMatch,floatMatch), t ).groups()

        info[ (int(matches[1]), int(matches[2]) ) ] = ( float(matches[8]), float(matches[9]), float(matches[10]) )

    return info

def readSModelXsection( filename ):
    import re
    f = open( filename )

    xSections = {}

    for line in f.readlines():
        if line.startswith("#"): continue
        m, xsec, uncert = line.split(" ")
        xSections[int(m)] = (float(xsec), float(uncert) )
    f.close()

    return xSections

def getCrossSectionDict( filename, signalPoints ):
    # file name determines the config file
    # signalpoints (list) determines the keys
    outDict = {}
    if "T5" in filename:
        xsecs = readSModelXsection( crossSectionFileT5 )
        for signalPoint in signalPoints:
            pass
    if "Bino" in filename or "Wino" in filename:
        # GGM
        crossSectionFile = crossSectionFileBino if "Bino" in filename else crossSectionFileWino
        xsecs = readSignalXSection( crossSectionFile )
        for point in signalPoints:
            match = re.match( "Sq(\d+)_Gl(\d+)", point )
            if match:
                mSq, mGl = match.groups()
                outDict[ point ] = xsecs[ int(mSq), int(mGl ) ][0]
    return outDict


def getSignalPoints( filename ):
    f = ROOT.TFile( filename )
    outList = []
    for key in f.GetListOfKeys():
        name = key.GetName()
        if name.startswith("pdfTree"):
            outList.append( name.replace("pdfTree","") )
    return outList

def getPdfSets():
    pdfNames = [
        "CT10",
        "CT10nnlo_as_0117",
        "CT10nnlo_as_0119",
        "MSTW2008nnlo68cl",
        "MSTW2008nnlo68cl_asmz-68cl",
        "MSTW2008nnlo68cl_asmz+68cl",
        "NNPDF23_nnlo_as_0116",
        "NNPDF23_nnlo_as_0117",
        "NNPDF23_nnlo_as_0118",
        "NNPDF23_nnlo_as_0119",
        "NNPDF23_nnlo_as_0120",
        "NNPDF23_nnlo_as_0121",
        "NNPDF23_nnlo_as_0122"
    ]

    out = {}
    for name in pdfNames:
        out[name] = lhapdf.mkPDFs( name )
    return out

def getWeightedHistograms( filename, treename, nExp, pdfs ):
    tree = readTree( filename, treename )
    nGen = tree.GetEntries()
    originalPdf = lhapdf.mkPDF("cteq6l1")

    out = {}
    binning = [ 100, 120, 160, 200, 270, 350 ] # overflow ist last bin
    import array
    for name, pdfSet in pdfs.iteritems():
        for i in range(len(pdfSet)):
            out[name+"_%s"%i] = ROOT.TH1F( name+"_%s"%i, "", len(binning)-1, array.array('d', binning) )

    for event in tree:
        originalWeight = originalPdf.xfxQ( event.id1, event.x1, event.scale ) * originalPdf.xfxQ( event.id2, event.x2, event.scale )
        for name, pdfSet in pdfs.iteritems():
            for i,pdf in enumerate(pdfSet):
                weight = pdf.xfxQ( event.id1, event.x1, event.scale ) * pdf.xfxQ( event.id2, event.x2, event.scale ) / originalWeight
                out[name+"_%s"%i].Fill( event.met, event.weight*weight*nExp/nGen )
    return out

def mergeNNPDF( histos ):
    # create output histos
    h_mean = histos.itervalues().next().Clone()
    h_up = h_mean.Clone()
    h_down = h_mean.Clone()

    for h in h_mean, h_up, h_down:
        h.Reset("ICESM")

    scales = { # gaus
        116: 4,
        117: 25,
        118: 71,
        119: 100,
        120: 71,
        121: 25,
        122: 4
    }

    # To make the second loop faster
    orderedValues = [[ None for i in range(100)] for j in range(len(scales)) ]

    for bin in range( h_mean.GetNbinsX()+2 ):
        f_sum_events = 0
        f_sum_replica_weights = 0
        for iAlphaS, alphaS in enumerate(range( 116, 123 )):
            sum_events_nnpdf =0
            for ij, j in enumerate(range( 1, 101 )):
                value = histos[ "NNPDF23_nnlo_as_0%s_%s"%(alphaS,j) ].GetBinContent( bin )
                orderedValues[iAlphaS][ij] = value
                sum_events_nnpdf += scales[alphaS] * value
            f_sum_replica_weights += (scales[alphaS]*100)
            f_sum_events += sum_events_nnpdf
        f_nnpdf_mean = f_sum_events / ( f_sum_replica_weights )

        f_sum_diffsq_events = 0
        for iAlphaS, alphaS in enumerate(range( 116, 123 )):
            sum_diffsq_events_nnpdf =0
            for j in orderedValues[iAlphaS]:
                sum_diffsq_events_nnpdf += scales[alphaS] * ( j - f_nnpdf_mean )**2
            f_sum_diffsq_events += sum_diffsq_events_nnpdf
        f_sum_diffsq_events = math.sqrt( f_sum_diffsq_events / f_sum_replica_weights )

        h_mean.SetBinContent( bin, f_nnpdf_mean )
        h_up.SetBinContent( bin, f_nnpdf_mean + f_sum_diffsq_events )
        h_down.SetBinContent( bin, f_nnpdf_mean - f_sum_diffsq_events )
        for h in h_mean, h_up, h_down:
            h.SetBinError( bin, histos["NNPDF23_nnlo_as_0119_0"].GetBinError(bin) )

    return h_mean, h_up, h_down



def mergePDF( histos, nName, uName, dName, norm_pdf=1., norm_as_plus=1., norm_as_minus=1., as_plus_member=[1,0], as_minus_member=[2,0] ):
    # create output histos
    h_mean = histos.itervalues().next().Clone()
    h_up = h_mean.Clone()
    h_down = h_mean.Clone()

    for h in h_mean, h_up, h_down:
        h.Reset("ICESM")

    for bin in range( h_mean.GetNbinsX()+2 ):
        mean = histos[nName+"_0"].GetBinContent(bin)
        pdf_plus = 0
        pdf_minus = 0
        for i in range( 1, 1000, 2 ):
            try:
                up   = histos["%s_%s"%(nName,i)].GetBinContent(bin)
                down = histos["%s_%s"%(nName,i+1)].GetBinContent(bin)

                pdf_plus += ( max([up-mean, down-mean,0]) )**2
                pdf_minus += ( max([mean-up, mean-down,0]) )**2
            except:
                break

        pdf_plus = math.sqrt( pdf_plus ) / norm_pdf
        pdf_minus = math.sqrt( pdf_minus ) / norm_pdf

        as_plus = 1./norm_as_plus*(histos[uName+"_0"].GetBinContent(bin)-mean)
        as_minus = 1./norm_as_minus*(histos[dName+"_0"].GetBinContent(bin)-mean)

        statError = histos[nName+"_0"].GetBinError( bin )
        tot_plus = math.sqrt( pow( pdf_plus, 2 ) + pow( as_plus, 2 ) )
        tot_minus = math.sqrt( pow( pdf_minus, 2 ) + pow( as_minus, 2 ) )

        h_mean.SetBinContent( bin, mean )
        h_up  .SetBinContent( bin, max( 0, mean + tot_plus  ) )
        h_down.SetBinContent( bin, max( 0, mean - tot_minus ) )

        for h in h_mean, h_up, h_down:
            h.SetBinError( bin, statError )

    return h_mean, h_up, h_down



def writeToFile( histos ):
    f = ROOT.TFile( "pdfhistos_singlepoint.root", "RECREATE" )
    f.cd()
    for name, hist in histos.iteritems():
        hist.Write()
    f.Close()


def mergePdfHistos( histos ):

    hs_NNPDF = mergeNNPDF( histos )

    hs_CT = mergePDF(
        histos,
        "CT10",
        "CT10nnlo_as_0117",
        "CT10nnlo_as_0119",
        1.64485362695147308,
        8.23893630338557559e-01,
        8.23893630338557559e-01
    )
    hs_MSTW = mergePDF(
        histos,
        "MSTW2008nnlo68cl",
        "MSTW2008nnlo68cl_asmz+68cl",
        "MSTW2008nnlo68cl_asmz-68cl",
        norm_as_minus=1.25356543847045088
    )

    for h in hs_CT:
        h.SetLineColor(2)
    for h in hs_MSTW:
        h.SetLineColor(3)
    hs_NNPDF[0].SetMaximum(900)
    hs_NNPDF[0].Draw("hist, e") # blue
    hs_CT[0].Draw("same hist, e") # red
    hs_MSTW[0].Draw("same hist e") # green
    for h in hs_NNPDF[1:]: h.Draw("same hist ")
    for h in hs_MSTW[1:]: h.Draw("same hist ")
    for h in hs_CT[1:]: h.Draw("same hist ")

    ROOT.gPad.SaveAs("test.pdf")

    return [ hs_NNPDF, hs_CT, hs_MSTW ]

def totalMinMaxHisto( histos, errors=True ):
    mergedList = []
    for l in histos:
        mergedList.extend( l )

    h = histos[0][0].Clone()
    h.Reset("ICESM")

    for bin in range( h.GetNbinsX()+2 ):
        mini = min( [ x.GetBinContent(bin) for x in mergedList ] )
        if errors:
            miniError = min( [ x.GetBinContent(bin)-x.GetBinError(bin) for x in mergedList ] )
            mini = min( [ miniError, mini ] )

        maxi = max( [ x.GetBinContent(bin) for x in mergedList ] )
        if errors:
            maxiError = max( [ x.GetBinContent(bin)+x.GetBinError(bin) for x in mergedList ] )
            maxi = max( [ maxiError, maxi ] )

        h.SetBinContent( bin, (mini+maxi)/2 )
        h.SetBinError( bin, (maxi-mini)/2 )
        if maxi+mini>0:
#            print (maxi-mini)/(maxi+mini)
            print math.sqrt( ((maxi-mini)/(maxi+mini))**2 - 0.196208658373**2 )
        else:
            print 0, "bin", bin
    return h


if __name__ == "__main__":
    filename = "GGM_Bino_V04.01_tree.root"
    signalPoints = getSignalPoints( filename )
    xSections = getCrossSectionDict( filename, signalPoints )
    nExpected = {}
    for sPoint, xsec in xSections.iteritems():
        nExpected[sPoint] = xsec*lumiData

    #pdfs = getPdfSets()

    for signalPoint in signalPoints:
        #histos = getWeightedHistograms( filename, "pdfTree"+signalPoint, nExpected[signalPoint], pdfs )
        #pickle.dump( histos, open( "histogram.pickle", "wb" ) )
        histos = pickle.load( open( "histogram.pickle", "rb" ) )
        binning = [ -10,10000 ] # overflow ist last bin
        import array
        #for name, h in histos.iteritems():
        #    histos[name] = h.Rebin(len(binning)-1, name, array.array("d",binning  ) )
        #writeToFile( histos )

        # merge histos
        mergedHistos = mergePdfHistos( histos )
        totalMinMax = totalMinMaxHisto( mergedHistos, errors=True )
        totalMinMax.Draw("hist e")
        ROOT.gPad.SaveAs("test2.pdf")


        import sys
        sys.exit()




