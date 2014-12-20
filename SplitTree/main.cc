#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TPRegexp.h"

std::string getModelFromGridParamStr( const std::vector<std::string>& gridParamStr ) {
    if( gridParamStr.size() < 2 ) return "";
    std::string line = gridParamStr.at(1);
    line = line.substr( 8 ); //all from T5.. on
    int pos = line.find( " " ); // find first space
    line = line.substr( 0, pos );
    return line;
}

std::string random_string( const int len ) {
    /* Creates a random string of length 'len'. */
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    std::string s( len, 'x' );
    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return s;
}

int main( int argc, char**argv ) {

    // set new random seed for random names
    srand(time(0));

    // Set error level
    //gErrorIgnoreLevel = kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;

    if( argc != 2 ) {
        std::cout << "Usage: " << argv[0] << " inputFilename.root" << std::endl;
        return 1;
    }
    std::cout << "Inputfilename: " << argv[1] << std::endl;

    // Find suffix of input file (will be appended to output filename)
    TPRegexp expFilename( ".*susyEvents(.*).root" );
    TObjArray *arr = expFilename.MatchS( argv[1] );
    std::string inputFileSuffix("");
    if( arr->GetLast() > 0 ){
        inputFileSuffix = (std::string)(((TObjString *)arr->At(1))->GetString());
    }

    // Read input tree
    std::string dCachePrefix( "dcap://dcache-cms-dcap.desy.de" );
    TChain oldTree( "susyTree" );
    if( !oldTree.Add( (dCachePrefix+argv[1]).c_str() ) ) {
        std::cout << "File could not be opened or tree could not be found" << std::endl;
        return 0;
    }

    // object to read
    std::vector<std::string>* gridParamStr = 0;
    oldTree.SetBranchAddress( "gridParamStr", &gridParamStr );

    std::string currentModel = "";
    std::string oldModel = "";
    TFile* newFile;
    TTree* newTree;

    // Start the event loop
    for( int i=0; i<oldTree.GetEntries(); i++) {
        oldTree.GetEntry(i);
        if( gridParamStr->size() != 2 ) {
            std::cout << "ERROR: gridParamStr must be 2" << std::endl;
            continue;
        }
        currentModel = getModelFromGridParamStr( *gridParamStr );
        if( currentModel.empty() ) std::cout << "ERROR: Something went wrong with extracting grid point" << std::endl;

        if( currentModel != oldModel ) {

            if ( i ) { // Write oldModel to file, but not for first entry
                std::cout << "Writing " << newTree->GetEntries() << " events to file " << newFile->GetName() << std::endl;
                newFile->Write();
                std::cout << "New model: " << currentModel << std::endl;
            } else {
                std::cout << "Beginning with point: " << currentModel << std::endl;
            }

            newFile = new TFile( ("/nfs/dust/cms/user/kiesel/tmp/"+currentModel+inputFileSuffix+"_"+random_string(9)+".root").c_str(), "NEW" ); // or maybe substr from in put file
            newTree = oldTree.CloneTree(0);
            oldModel = currentModel;
        }
        newTree->Fill();
    }
    newFile->Write();
    delete gridParamStr;
    delete newTree;
    delete newFile;
}

