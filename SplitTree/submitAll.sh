#!/bin/bash

datasetsGG=(
T5gg_800_V04
T5gg_1000_V04
T5gg_1200_V04
T5gg_1350_V04
T5gg_1550_V04
)
datasetsGW=(
T5wg_400_V04
T5wg_800_V04
T5wg_600_V04
T5wg_1000_V04
T5wg_1200_V04
T5wg_1350_V04
T5wg_1550_V04
)

submit=false

for dataset in "${datasetsGW[@]}"; do
    for file in $(srmls srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=/pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/$dataset/|grep root |awk '{print $2}'); do
        jobName1=$(echo $file|cut -d '/' -f10)
        jobName2=$(echo $file|cut -d '_' -f4)
        if [[ $file == *susyEvents_57_1_ama.root ]]; then
            submit=true
        fi
        #qsub -N ${jobName1}-${jobName2} -b y -j y -l os=sld6 -l h_vmem=1000M -l h_rt=1:00:00 -l site=hh  /afs/desy.de/user/k/kiesel/scratch/singlePhoton/SplitTree/submit $file
        if [ "$submit" == true ]; then
            /afs/desy.de/user/k/kiesel/scratch/singlePhoton/SplitTree/submit $file
        fi
    done
    exit # for testing reasons only one set for now
done
