#!/usr/bin/env python
import ROOT
import re
import lhapdf

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
        out[name] = lhapdf.mkPDF( name, 0 )
    return out

def getWeightedHistograms( filename, treename, nExp, pdfs ):
    tree = readTree( filename, treename )
    nGen = tree.GetEntries()
    originalPdf = lhapdf.mkPDF("cteq6l1")

    out = {}
    binning = [ 100, 120, 160, 200, 270, 350, 2000 ]
    import array
    for name in pdfs:
        out[name] = ROOT.TH1F( name, "", len(binning)-1, array.array('d', binning) )

    for event in tree:
        originalWeight = originalPdf.xfxQ( event.id1, event.x1, event.scale ) * originalPdf.xfxQ( event.id2, event.x2, event.scale )
        for name, pdf in pdfs.iteritems():
            weight = pdf.xfxQ( event.id1, event.x1, event.scale ) * pdf.xfxQ( event.id2, event.x2, event.scale ) / originalWeight
            out[name].Fill( event.met, event.weight*weight*nExp/nGen )
    return out

def mergeNNPDF( histos ):
    # Just for make life simpler: copy the important histograms and rename the keys
    # The keys are 1/alpha_s (all around 119)
    h = {}
    for name, histo in histos.iteritems():
        if name.startswith( "NNPDF23_nnlo_as_0" ):
            newName = int(name.replace( "NNPDF23_nnlo_as_0", "" ))
            h[newName] = histo

    # create output histos
    h_mean = h[119].Clone()
    h_up__ = h[119].Clone()
    h_down = h[119].Clone()

    for h in h_mean, h_up__, h_down:
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
    for bin in range( h_mean.GetNbinsX()+2 ):
        bin_mean = h[119]->GetBinContent(bin)

    """for bin in range( h_mean.GetNbinsX()+2 ):
        bin_mean = h[119]->GetBinContent(bin)
        pdf_plus  = 0
        pdf_minus = 0
        for index in range( 0.5*(len(h)-1) ):
            up   = h[119+index]->GetBinContent(bin)
            down = h[119-index]->GetBinContent(bin)
            pdf_plus  += pow( max([up-bin_mean, down-bin_mean, 0]), 2 )
            pdf_minus += pow( max([bin_mean-up, bin_mean-down, 0]), 2 )

        ??? 16?
    """


def mergePDF( nHisto, uHisto, dHisto, norm_pdf=1., norm_as_plus=1., norm_as_minus=1., as_plus_member=[1,0], as_minus_member=[2,0] ):

    # create output histos
    h_up   = nHisto.Clone()
    h_down = nHisto.Clone()

    for h in h_up, h_down:
        h.Reset("ICESM")

    for bin in range( nHisto.GetNbinsX()+2 ):
        mean = nHisto.GetBinContent(bin)
        up   = uHisto.GetBinContent(bin)
        down = dHisto.GetBinContent(bin)

        pdf_plus = max([up-mean, down-mean,0]) / norm_pdf
        pdf_minus = min([mean-up, mean-down,0] / norm_pdf

        as_plus = 1./norm_as_plus*(up-mean)
        as_minus = 1./norm_as_minus*(down-mean)

        statError = nHisto.GetBinError( bin )

        tot_plus = sqrt( pow( pdf_plus, 2 ) + pow( as_plus, 2 ) )
        tot_minus = sqrt( pow( pdf_minus, 2 ) + pow( as_minus, 2 ) )

        h_up  .SetBinContent( bin, max( 0, mean + tot_plus  ) )
        h_down.SetBinContent( bin, max( 0, mean + tot_minus ) )

        h_up  .SetBinError( bin, statError )
        h_down.SetBinError( bin, statError )

    return nHisto, h_up, h_down








def mergePdfHistos( histos ):

    h_NNPDF = mergeNNPDF( histos )
    h_CT = mergePDF(
        histos["CT10"],
        histos["CT10nnlo_as_0117"],
        histos["CT10nnlo_as_0119"],
        1.64485362695147308,
        8.23893630338557559e-01,
        8.23893630338557559e-01
    )
    h_MSTW = mergePDF(
        histos["MSTW2008nnlo68cl"],
        histos["MSTW2008nnlo68cl_asmz+68cL"],
        histos["MSTW2008nnlo68cl_asmz-68cl"],
        norm_as_minus=1.25356543847045088
    )

    return histos


if __name__ == "__main__":
    filename = "GGM_Bino_V04.01_tree.root"
    signalPoints = getSignalPoints( filename )
    xSections = getCrossSectionDict( filename, signalPoints )
    nExpected = {}
    for sPoint, xsec in xSections.iteritems():
        nExpected[sPoint] = xsec*lumiData

    pdfs = getPdfSets()

    for signalPoint in signalPoints:
        histos = getWeightedHistograms( filename, "pdfTree"+signalPoint, nExpected[signalPoint], pdfs )

        # merge histos
        mergedHistos = mergePdfHistos( histos )


        import sys
        sys.exit()




