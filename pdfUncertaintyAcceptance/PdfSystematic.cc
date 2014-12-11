#include "LHAPDF/LHAPDF.h"

#include<TTree.h>
#include<TFile.h>
#include<TF1.h>

#include<iostream>
#include<vector>

using namespace std;

typedef vector<string> strings;

class PdfSystematicsAnalyzer {
    // This class is a clone of CMSSW/ElectroWeakAnalysis/Utilities/src/PdfWeightProducer.cc
    // For more information visit
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEWKUtilities#PDF_SYSTEMATICS_PDFWeightProduce

    public:
        // Normally, one object of this class is initialized
        PdfSystematicsAnalyzer( const strings &pdfNames );
        ~PdfSystematicsAnalyzer();

        // Pdf analysis for this tree. All internal variables from
        // other trees are reset.
        void newTree( TTree &tree );

        // These two functions have to be called for each tree.
        void fillWeights();
        void calculate();

    private:
        // Input tree
        TTree* tree_;

        strings pdfNames_;
        vector< vector<LHAPDF::PDF*> > pdfSets_;
        unsigned int originalEvents_;
        unsigned int selectedEvents_;
        std::vector<int> pdfStart_;
        std::vector<double> weightedSelectedEvents_;
        std::vector<double> weighted2SelectedEvents_;
        std::vector<double> weightedEvents_;

        // Needed for merging nnpdf uncertainties
        float gausWeights_[4];
};


PdfSystematicsAnalyzer::PdfSystematicsAnalyzer( const strings &pdfNames ) :
    pdfNames_(pdfNames)
{
    // Fills the array gausWeights with 1, 0.61, 0.14, 0.01
    TF1 gaus("gaus", "gaus", 0, 4);
    gaus.SetParameters(1,0,1);
    float sumW = 0;
    for (unsigned int i=0; i<4; i++) {
        gausWeights_[i] = gaus.Eval( i );
        sumW += gausWeights_[i];
        if(i>0) sumW += gausWeights_[i]; // all entries except the center one will be used on each side (2 times)
    }

    // Normalize the array
    for (unsigned int i=0; i<4; i++) {
        gausWeights_[i] /= sumW;
    }

    for (unsigned int i=0; i<pdfNames_.size(); ++i) {
        LHAPDF::PDFSet set( pdfNames[i] );
        pdfSets_.push_back( set.mkPDFs() );
    }
}

void PdfSystematicsAnalyzer::newTree( TTree &tree )
{
    tree_ = &tree;
    originalEvents_ = 0;
    selectedEvents_ = 0;
    pdfStart_.clear();
    weightedSelectedEvents_.clear();
    weighted2SelectedEvents_.clear();
    weightedEvents_.clear();

    for (unsigned int i=0; i<pdfNames_.size(); ++i) {
        pdfStart_.push_back(-1);
    }
}

/////////////////////////////////////////////////////////////////////////////////////
PdfSystematicsAnalyzer::~PdfSystematicsAnalyzer()
{
}

/////////////////////////////////////////////////////////////////////////////////////
void PdfSystematicsAnalyzer::fillWeights(){

    // Initialize tree
    float x1, x2, scale, weight, met, id1, id2;

    tree_->SetBranchAddress( "met", &met );
    tree_->SetBranchAddress( "x1", &x1 );
    tree_->SetBranchAddress( "x2", &x2 );
    tree_->SetBranchAddress( "scale", &scale );
    tree_->SetBranchAddress( "weight", &weight );
    tree_->SetBranchAddress( "id1", &id1 );
    tree_->SetBranchAddress( "id2", &id2 );

    originalEvents_ = tree_->GetEntries();

    for( unsigned int entry=0; entry< originalEvents_; ++entry ) {
        tree_->GetEntry( entry );

        // Here the cut is defined
        bool selectedEvent = met >= 100;

        if (selectedEvent) selectedEvents_++;

        for (int i=0; i<pdfNames_.size(); ++i ) {
            unsigned int nmembers = pdfSets_[i].size();

            // Set up arrays the first time weights are read
            if (pdfStart_[i]<0) {
                  pdfStart_[i] = weightedEvents_.size();
                  for (unsigned int j=0; j<nmembers; ++j) {
                        weightedEvents_.push_back(0.);
                        weightedSelectedEvents_.push_back(0.);
                        weighted2SelectedEvents_.push_back(0.);
                  }
            }

            // Fill the weights
            for (unsigned int j=0; j<nmembers; ++j) {
                float weight = pdfSets_[i][j]->xfxQ( (int)id1, x1, scale ) * pdfSets_[i][j]->xfxQ( (int)id2, x2, scale );
                weightedEvents_[pdfStart_[i]+j] += weight;
                if (selectedEvent) {
                    weightedSelectedEvents_[pdfStart_[i]+j] += weight;
                    weighted2SelectedEvents_[pdfStart_[i]+j] += weight*weight;
                }
            }

        } // pdfnames
    } // entries
}

/////////////////////////////////////////////////////////////////////////////////////
void PdfSystematicsAnalyzer::calculate(){
    int log = 0;

    if (originalEvents_==0) {
        cout << "NO EVENTS => NO RESULTS" << endl;
        return;
    }
    if (selectedEvents_==0) {
        cout << "NO SELECTED EVENTS => NO RESULTS" << endl;
        return;
    }

    if(log>4) cout << "\n>>>> Begin of PDF weight systematics summary >>>>" << endl;
    if(log>4) cout << "Total number of analyzed data: " << originalEvents_ << " [events]" << endl;
    double originalAcceptance = double(selectedEvents_)/originalEvents_;
    if(log>4) cout << "Total number of selected data: " << selectedEvents_ << " [events], corresponding to acceptance: [" << originalAcceptance*100 << " +- " << 100*sqrt( originalAcceptance*(1.-originalAcceptance)/originalEvents_) << "] %" << endl;

    if(log>4) cout << "\n>>>>> PDF UNCERTAINTIES ON RATE >>>>>>" << endl;
    for (unsigned int i=0; i<pdfNames_.size(); ++i) {
        bool nnpdfFlag = ( pdfNames_[i].find("NNPDF") != string::npos );
        unsigned int nmembers = weightedSelectedEvents_.size()-pdfStart_[i];
        if (i<pdfNames_.size()-1) nmembers = pdfStart_[i+1] - pdfStart_[i];
        unsigned int npairs = (nmembers-1)/2;
        if(log>4) cout << "RATE Results for PDF set " << pdfNames_[i] << " ---->" << endl;

        double events_central = weightedSelectedEvents_[pdfStart_[i]];
        if(log>4) cout << "\tEstimate for central PDF member: " << int(events_central) << " [events]" << endl;
        double events2_central = weighted2SelectedEvents_[pdfStart_[i]];
        if(log>4) cout << "\ti.e. [" << std::setprecision(4) << 100*(events_central-selectedEvents_)/selectedEvents_ << " +- " <<
            100*sqrt(events2_central-events_central+selectedEvents_*(1-originalAcceptance))/selectedEvents_
        << "] % relative variation with respect to original PDF" << endl;

        if (npairs>0) {
            if(log>4) cout << "\tNumber of eigenvectors for uncertainty estimation: " << npairs << endl;
            double wplus = 0.;
            double wminus = 0.;
            unsigned int nplus = 0;
            unsigned int nminus = 0;
            for (unsigned int j=0; j<npairs; ++j) {
                double wa = weightedSelectedEvents_[pdfStart_[i]+2*j+1]/events_central-1.;
                double wb = weightedSelectedEvents_[pdfStart_[i]+2*j+2]/events_central-1.;
                if (nnpdfFlag) {
                    if (wa>0.) {
                        wplus += wa*wa;
                        nplus++;
                    } else {
                        wminus += wa*wa;
                        nminus++;
                    }
                    if (wb>0.) {
                        wplus += wb*wb;
                        nplus++;
                    } else {
                        wminus += wb*wb;
                        nminus++;
                    }
                } else {
                    if (wa>wb) {
                        if (wa<0.) wa = 0.;
                        if (wb>0.) wb = 0.;
                        wplus += wa*wa;
                        wminus += wb*wb;
                    } else {
                        if (wb<0.) wb = 0.;
                        if (wa>0.) wa = 0.;
                        wplus += wb*wb;
                        wminus += wa*wa;
                    }
                }
            }
            if (wplus>0) wplus = sqrt(wplus);
            if (wminus>0) wminus = sqrt(wminus);
            if (nnpdfFlag) {
               if (nplus>0) wplus /= sqrt(nplus);
               if (nminus>0) wminus /= sqrt(nminus);
            }
            if(log>4) cout << "\tRelative uncertainty with respect to central member: +" << std::setprecision(4) << 100.*wplus << " / -" << std::setprecision(4) << 100.*wminus << " [%]" << endl;
        } else {
            if(log>4) cout << "\tNO eigenvectors for uncertainty estimation" << endl;
        }
    }

    // Variables needed for combining alpha_s uncertainties
    float ctUp(0), ctDn(0), ctAccAlphaS_Up(0), ctAccAlphaS_Dn(0);
    float mstwUp(0), mstwDn(0), mstwAccAlphaS_Up(0), mstwAccAlphaS_Dn(0);
    float nnpdfUp(0), nnpdfDn(0), nnpdfAccAlphaS_116(0), nnpdfAccAlphaS_117(0), nnpdfAccAlphaS_118(0), nnpdfAccAlphaS_119(0), nnpdfAccAlphaS_120(0), nnpdfAccAlphaS_121(0), nnpdfAccAlphaS_122(0);

    if(log>4) cout << "\n>>>>> PDF UNCERTAINTIES ON ACCEPTANCE >>>>>>" << endl;
    for (unsigned int i=0; i<pdfNames_.size(); ++i) {
        bool nnpdfFlag = ( pdfNames_[i].find("NNPDF") != string::npos );
        unsigned int nmembers = weightedEvents_.size()-pdfStart_[i];
        if (i<pdfNames_.size()-1) nmembers = pdfStart_[i+1] - pdfStart_[i];
        unsigned int npairs = (nmembers-1)/2;
        if(log>4) cout << "ACCEPTANCE Results for PDF set " << pdfNames_[i] << " ---->" << endl;

        double acc_central = 0.;
        double acc2_central = 0.;
        if (weightedEvents_[pdfStart_[i]]>0) {
            acc_central = weightedSelectedEvents_[pdfStart_[i]]/weightedEvents_[pdfStart_[i]];
            acc2_central = weighted2SelectedEvents_[pdfStart_[i]]/weightedEvents_[pdfStart_[i]];
        }

        // Fill acceptances
        if (pdfNames_[i] == "CT10nnlo_as_0117" ) {
            ctAccAlphaS_Up = acc_central;
        } else if (pdfNames_[i] == "CT10nnlo_as_0119" ) {
            ctAccAlphaS_Dn = acc_central;
        } else if (pdfNames_[i] == "MSTW2008nnlo68cl_asmz+68cl" ) {
            mstwAccAlphaS_Up = acc_central;
        } else if (pdfNames_[i] == "MSTW2008nnlo68cl_asmz-68cl" ) {
            mstwAccAlphaS_Dn = acc_central;
        } else if (pdfNames_[i] == "NNPDF23_nnlo_as_0116" ) {
            nnpdfAccAlphaS_116 = acc_central;
        } else if (pdfNames_[i] == "NNPDF23_nnlo_as_0117" ) {
            nnpdfAccAlphaS_117 = acc_central;
        } else if (pdfNames_[i] == "NNPDF23_nnlo_as_0118" ) {
            nnpdfAccAlphaS_118 = acc_central;
        } else if (pdfNames_[i] == "NNPDF23_nnlo_as_0119" ) {
            nnpdfAccAlphaS_119 = acc_central;
        } else if (pdfNames_[i] == "NNPDF23_nnlo_as_0120" ) {
            nnpdfAccAlphaS_120 = acc_central;
        } else if (pdfNames_[i] == "NNPDF23_nnlo_as_0121" ) {
            nnpdfAccAlphaS_121 = acc_central;
        } else if (pdfNames_[i] == "NNPDF23_nnlo_as_0122" ) {
            nnpdfAccAlphaS_122 = acc_central;
        }

        double waverage = weightedEvents_[pdfStart_[i]]/originalEvents_;
        if(log>4) cout << "\tEstimate for central PDF member acceptance: [" << acc_central*100 << " +- " <<
        100*sqrt((acc2_central/waverage-acc_central*acc_central)/originalEvents_)
        << "] %" << endl;
        double xi = acc_central-originalAcceptance;
        double deltaxi = (acc2_central-(originalAcceptance+2*xi+xi*xi))/originalEvents_;
        if (deltaxi>0) deltaxi = sqrt(deltaxi); //else deltaxi = 0.;
        if(log>4) cout << "\ti.e. [" << std::setprecision(4) << 100*xi/originalAcceptance << " +- " << std::setprecision(4) << 100*deltaxi/originalAcceptance << "] % relative variation with respect to the original PDF" << endl;

        if (npairs>0) {
            if(log>4) cout << "\tNumber of eigenvectors for uncertainty estimation: " << npairs << endl;
          double wplus = 0.;
          double wminus = 0.;
          unsigned int nplus = 0;
          unsigned int nminus = 0;
          for (unsigned int j=0; j<npairs; ++j) {
            double wa = 0.;
            if (weightedEvents_[pdfStart_[i]+2*j+1]>0) wa = (weightedSelectedEvents_[pdfStart_[i]+2*j+1]/weightedEvents_[pdfStart_[i]+2*j+1])/acc_central-1.;
            double wb = 0.;
            if (weightedEvents_[pdfStart_[i]+2*j+2]>0) wb = (weightedSelectedEvents_[pdfStart_[i]+2*j+2]/weightedEvents_[pdfStart_[i]+2*j+2])/acc_central-1.;
            if (nnpdfFlag) {
                if (wa>0.) {
                    wplus += wa*wa;
                    nplus++;
                } else {
                    wminus += wa*wa;
                    nminus++;
                }
                if (wb>0.) {
                    wplus += wb*wb;
                    nplus++;
                } else {
                    wminus += wb*wb;
                    nminus++;
                }
            } else {
                if (wa>wb) {
                    if (wa<0.) wa = 0.;
                    if (wb>0.) wb = 0.;
                    wplus += wa*wa;
                    wminus += wb*wb;
                } else {
                    if (wb<0.) wb = 0.;
                    if (wa>0.) wa = 0.;
                    wplus += wb*wb;
                    wminus += wa*wa;
                }
            }
          }
          if (wplus>0) wplus = sqrt(wplus);
          if (wminus>0) wminus = sqrt(wminus);
          if (nnpdfFlag) {
            if (nplus>0) wplus /= sqrt(nplus);
            if (nminus>0) wminus /= sqrt(nminus);
          }
          if(log>4) cout << "\tRelative uncertainty with respect to central member: +" << std::setprecision(4) << 100.*wplus << " / -" << std::setprecision(4) << 100.*wminus << " [%]" << endl;

          // Fill variables (2)
          if( pdfNames_[i] == "CT10" ) {
              ctUp = wplus;
              ctDn = wminus;
          } else if( pdfNames_[i] == "MSTW2008nnlo68cl" ) {
              mstwUp = wplus;
              mstwDn = wminus;
          } else if( pdfNames_[i] == "NNPDF23_nnlo_as_0119" ) {
              nnpdfUp = wplus;
              nnpdfDn = wminus;
          }

        } else {
            if(log>4) cout << "\tNO eigenvectors for uncertainty estimation" << endl;
        }
    }
    if(log>4) cout << ">>>> End of PDF weight systematics summary >>>>" << endl;

    // add quadratically uncertanity of pdf and alpha_S
    // for CT10, the uncertainty is 90% cl, so it has to be reduced by a factor of 1.645
    float ctTotalUp = (ctUp-ctDn)/2 + sqrt( pow((ctUp+ctDn)/2,2) + pow((ctAccAlphaS_Up-ctAccAlphaS_Dn)/2,2) )/1.645;
    float ctTotalDn =-(ctUp-ctDn)/2 + sqrt( pow((ctUp+ctDn)/2,2) + pow((ctAccAlphaS_Up-ctAccAlphaS_Dn)/2,2) )/1.645;

    float mstwTotalUp = (mstwUp-mstwDn)/2 + sqrt( pow((mstwUp+mstwDn)/2,2) + pow((mstwAccAlphaS_Up-mstwAccAlphaS_Dn)/2,2) );
    float mstwTotalDn =-(mstwUp-mstwDn)/2 + sqrt( pow((mstwUp+mstwDn)/2,2) + pow((mstwAccAlphaS_Up-mstwAccAlphaS_Dn)/2,2) );

    float nnpdfAlphaS_UncertSquare =
        gausWeights_[3] * pow( nnpdfAccAlphaS_119 - nnpdfAccAlphaS_116, 2 ) +
        gausWeights_[2] * pow( nnpdfAccAlphaS_119 - nnpdfAccAlphaS_117, 2 ) +
        gausWeights_[1] * pow( nnpdfAccAlphaS_119 - nnpdfAccAlphaS_118, 2 ) +
        gausWeights_[1] * pow( nnpdfAccAlphaS_119 - nnpdfAccAlphaS_120, 2 ) +
        gausWeights_[2] * pow( nnpdfAccAlphaS_119 - nnpdfAccAlphaS_121, 2 ) +
        gausWeights_[3] * pow( nnpdfAccAlphaS_119 - nnpdfAccAlphaS_122, 2 );
    float nnpdfTotalUp = (nnpdfUp-nnpdfDn)/2 + sqrt( pow((nnpdfUp+nnpdfDn)/2,2) + nnpdfAlphaS_UncertSquare );
    float nnpdfTotalDn =-(nnpdfUp-nnpdfDn)/2 + sqrt( pow((nnpdfUp+nnpdfDn)/2,2) + nnpdfAlphaS_UncertSquare );

    float totalUp = max( max( ctTotalUp, mstwTotalUp ), nnpdfTotalUp );
    float totalDn = max( max( ctTotalDn, mstwTotalDn ), nnpdfTotalDn );

    cout << tree_->GetName() << " " << (totalUp+totalDn)/2 << endl;


}

int main (int argc, char** argv) {
    if( argc < 2 ) {
        cerr << "Usage: " <<argv[0] << " filename1.root [filename2.root]" << endl;
        return 1;
    }

    strings pdfNames;
    pdfNames.push_back( "CT10" );
    pdfNames.push_back( "CT10nnlo_as_0117" );
    pdfNames.push_back( "CT10nnlo_as_0119" );

    pdfNames.push_back( "MSTW2008nnlo68cl" );
    pdfNames.push_back( "MSTW2008nnlo68cl_asmz+68cl" );
    pdfNames.push_back( "MSTW2008nnlo68cl_asmz-68cl" );

    pdfNames.push_back( "NNPDF23_nnlo_as_0116" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0117" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0118" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0119" ); // main pdf
    pdfNames.push_back( "NNPDF23_nnlo_as_0120" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0121" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0122" );

   PdfSystematicsAnalyzer analyzer( pdfNames );

    for( int i=1; i<argc; i++ ) {
        TFile file( argv[i] );
        TIter next( file.GetListOfKeys() );
        TObject* obj;
        while (( obj = next() ) ) {
            string objName = obj->GetName();
            if( objName.find("pdfTree") != std::string::npos ) {
                analyzer.newTree( *((TTree*)file.Get( objName.c_str() )));
                analyzer.fillWeights();
                analyzer.calculate();
            }
        }
        file.Close();
    }
}



















