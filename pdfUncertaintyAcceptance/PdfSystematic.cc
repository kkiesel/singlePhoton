#include "LHAPDF/LHAPDF.h"

#include<TTree.h>
#include<TFile.h>

#include<iostream>
#include<vector>

using namespace std;
typedef vector<string> strings;

class PdfSystematicsAnalyzer {
    public:
        PdfSystematicsAnalyzer( const strings &pdfNames );
        ~PdfSystematicsAnalyzer();
        void newTree( TTree &tree );
        void fillWeights();
        void calculate();
    private:
        TTree* tree_;
        strings pdfNames_;
        vector< vector<LHAPDF::PDF*> > pdfSets_;
        unsigned int originalEvents_;
        unsigned int selectedEvents_;
        std::vector<int> pdfStart_;
        std::vector<double> weightedSelectedEvents_;
        std::vector<double> weighted2SelectedEvents_;
        std::vector<double> weightedEvents_;
};


PdfSystematicsAnalyzer::PdfSystematicsAnalyzer( const strings &pdfNames ) :
    pdfNames_(pdfNames)
{
    for (unsigned int i=0; i<pdfNames_.size(); ++i) {
        LHAPDF::PDFSet set( pdfNames[i] );
        pdfSets_.push_back( set.mkPDFs() );
    }
};

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
    delete tree_;
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

    float totalUp = 0;
    float totalDown = 0;

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
            if( totalUp < wplus ) totalUp = wplus;
            if( totalDown < wminus ) totalDown = wminus;
        } else {
            if(log>4) cout << "\tNO eigenvectors for uncertainty estimation" << endl;
        }
    }
    if(log>4) cout << ">>>> End of PDF weight systematics summary >>>>" << endl;


    cout << tree_->GetName() << "\t" << (totalUp+totalDown)*50 << endl;

}

int main (int argc, char** argv) {
    if( argc < 2 ) {
        cerr << "Usage: " <<argv[0] << " filename1.root [filename2.root]" << endl;
        return 1;
    }

    strings pdfNames;
    pdfNames.push_back( "CT10" );
//    pdfNames.push_back( "CT10nnlo_as_0117" );
//    pdfNames.push_back( "CT10nnlo_as_0119" );

    pdfNames.push_back( "MSTW2008nnlo68cl" );
//    pdfNames.push_back( "MSTW2008nnlo68cl_asmz+68cl" );
//    pdfNames.push_back( "MSTW2008nnlo68cl_asmz-68cl" );

    pdfNames.push_back( "NNPDF23_nnlo_as_0119" );
/*    pdfNames.push_back( "NNPDF23_nnlo_as_0116" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0117" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0118" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0116" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0120" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0121" );
    pdfNames.push_back( "NNPDF23_nnlo_as_0122" );
*/
   PdfSystematicsAnalyzer analyzer( pdfNames );

    for( int i=1; i<argc; i++ ) {
        TFile file( argv[i] );
        TIter next( file.GetListOfKeys() );
        TObject* obj;
        while (( obj = next() ) ) {
            string objName = obj->GetName();
            if( objName.find("pdfTree") != std::string::npos ) {
                TTree* tree = (TTree*)file.Get( objName.c_str() );
                analyzer.newTree( *tree );
                analyzer.fillWeights();
                analyzer.calculate();
            }
        }
    }
}



















