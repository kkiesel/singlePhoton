#!/bin/bash -e
wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.1.0.tar.gz --no-check-certificate
tar -xf LHAPDF-6.1.0.tar.gz

wget http://yaml-cpp.googlecode.com/files/yaml-cpp-0.3.0.tar.gz -O- | tar xz
cd yaml-cpp
cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PWD/../local
make -j2 && make install
cd ..

cd LHAPDF-6.1.0

./configure --with-boost=/cvmfs/cms.cern.ch/slc6_amd64_gcc490/external/boost/1.51.0/ --prefix=$PWD/../local --with-yaml-cpp=$PWD/../local

# somehow the libpython2.7 is not found, a quick fix would be
# mkdir src/.libs
# ln -s /cvmfs/cms.cern.ch/slc6_amd64_gcc490/cms/cmssw/CMSSW_7_1_0_pre7/external/slc6_amd64_gcc490/lib/libpython2.7.so src/.libs/
make -j8
make install

cd ..
cd local/share/LHAPDF/
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0110.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0111.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0112.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0113.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0114.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0115.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0116.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0117.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0118.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0119.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0120.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0121.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0122.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0123.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0124.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0125.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0126.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0127.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0128.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0129.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0130.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10as.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/MSTW2008nnlo68cl.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/MSTW2008nnlo68cl_asmz+68cl.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/MSTW2008nnlo68cl_asmz-68cl.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/NNPDF23_nnlo_as_0116.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/NNPDF23_nnlo_as_0117.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/NNPDF23_nnlo_as_0118.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/NNPDF23_nnlo_as_0119.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/NNPDF23_nnlo_as_0120.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/NNPDF23_nnlo_as_0121.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/NNPDF23_nnlo_as_0122.tar.gz --no-check-certificate
wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/cteq6l1.tar.gz --no-check-certificate


for i in *.tar.gz
do
    tar -xf $i
done

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/local/lib/
export LHAPDF_BASE=$PWD/local/

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
