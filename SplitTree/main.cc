#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <vector>


std::string getModelFromGridParamStr( const std::vector<std::string>& gridParamStr ) {
    if( gridParamStr.size() < 2 ) return "";
    std::string line = gridParamStr.at(1);
    line = line.substr( 8 ); //all from T5.. on
    int pos = line.find( " " ); // find first space
    line = line.substr( 0, pos );
    return line;
}

int main( int argc, char** argv ) {
    if( argc < 2 ) {
        std::cout << "Usage: " << argv[0] << " inputFilename.root" << std::endl;
        return 1;
    }

    // creating input
    char* infilename = argv[1];

    // without .root
    std::string basename = infilename;
    basename = basename.substr( 0, basename.size()-5 );

    TChain oldTree( "susyTree" );
    oldTree.Add( infilename );

    // object to read
    std::vector<std::string>* gridParamStr = 0;
    oldTree.SetBranchAddress( "gridParamStr", &gridParamStr );

    std::string currentModel = "";
    std::string oldModel = "";
    TFile* newFile;
    TTree* newTree;

    for( int i=0; i<oldTree.GetEntries(); i++) {
        oldTree.GetEntry(i);
        if( gridParamStr->size() != 2 ) {
            std::cout << "ERROR: gridParamStr must be 2" << std::endl;
            continue;
        }
        currentModel = getModelFromGridParamStr( *gridParamStr );
        if( currentModel.empty() ) std::cout << "ERROR: Something went wrong with extracting grid point" << std::endl;

        if( currentModel != oldModel ) {
            std::cout << "New model: " << currentModel << std::endl;
            oldModel = currentModel;
            if( i ) newFile->Write(); // Write oldModel to file, but not for first entry
            newFile = new TFile( (basename+"_"+currentModel+".root").c_str(), "NEW" ); // or maybe substr from in put file
            newTree = oldTree.CloneTree(0);
        }
        newTree->Fill();
    }
    newFile->Write();
    delete gridParamStr;
    delete newTree;
    delete newFile;
}

