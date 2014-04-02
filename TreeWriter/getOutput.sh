#!/bin/bash
# Please call with version number to proceed as first option

datasets=(
PhotonHadA_V03
PhotonHadB_V03
PhotonHadC_V03
PhotonHadD_V03
GJets_200_400_V03
GJets_400_inf_V03
QCD_250_500_V03
QCD_500_1000_V03
QCD_1000_inf_V03
WJets_250_300_V03
WJets_300_400_V03
WJets_400_inf_V03
WJets_V02
TTJets_V03
WGamma_50_130_V03
WGamma_130_inf_V03
TTGamma_V03
ZGammaNuNu1_200_400_V03
ZGammaNuNu1_400_inf_V03
ZGammaNuNu2_200_400_V03
ZGammaNuNu2_400_inf_V03
ZGammaNuNu_V03
)
outputPath=/nfs/dust/cms/user/kiesel
version=$1
for dataset in "${datasets[@]}"; do

	abbr=$dataset.$version
	echo
	echo get out put for $abbr
	if ls $outputPath/${abbr}__*.root &> /dev/null; then
		hadd -f -k $outputPath/${abbr}_tree.root $outputPath/${abbr}__*.root
	fi

	# Clear information
	if ! qstat -u $USER -r|grep $abbr; then
		if ls ${abbr}__*.sh &> /dev/null; then
			tar -cfa $outputPath/logs/scripts_${abbr}.tar.bz2 ${abbr}__*.sh
		fi
		if ls ~/${abbr}__*.sh.o.* &> /dev/null; then
			tar -cfa $outputPath/logs/output_${abbr}.tar.bz2 ~/${abbr}__*.sh.o*
		fi
		rm -f ~/${abbr}__*.sh.o*
		rm -f ${abbr}__*.sh
		rm -f $outputPath/kiesel/${abbr}__*.root
	fi
done
