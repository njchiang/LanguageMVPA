#!/bin/bash
:<<doc
r to z transform was already run... so comparing those to their variance (sqrt N-3)
outputs:
n1000: randomising r output
avg_fdr: r --> z --> fdr
avg_fdrP: p(from analysis ) -->fdr
doc
sDir=/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/Maps/RSA
nDir=${sDir}/nMaps-9mm
mDir=/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/data/masks
FDRRATE=.05
n=1000 #randomise
zscale=1981 #number of items correlated
mkdir ${sDir}/zscored
mkdir ${sDir}/corrected
mkdir ${sDir}/output
cd ${sDir}/
for fileName in `ls nMaps-9mm/Searchlight_n17_nMap*`
do

maskName=`echo ${fileName} | cut -d'.' -f1 | sed 's/nMaps-9mm\/Searchlight_n17_nMap_//g'`

#scale up z
fslmaths ${fileName} -mul $zscale nMaps-9mm/${maskName}_denom #not sure if i actually need to multiply by number of matrices. makes sense though
fslmaths nMaps-9mm/${maskName}_denom -sub 3 nMaps-9mm/${maskName}_denom
fslmaths nMaps-9mm/${maskName}_denom -sqrt nMaps-9mm/${maskName}_denom

for modelName in Language Pictures LvP StimType Syntax AllSemantics Semantics LSemantics PSemantics LvPSemantics
do
echo "zscoring"
fslmaths concat/Searchlight_n17_rMap_${maskName}_${modelName} -mul nMaps-9mm/${maskName}_denom zscored/${maskName}_${modelName}
fslmaths zscored/${maskName}_${modelName} -ztop zscored/p_${maskName}_${modelName} #p values for each subject in case we want

echo "Averaging"
fslmaths concat/Searchlight_n17_rMap_${maskName}_${modelName} -Tmean output/avg_corr_${maskName}_${modelName}

echo "fdr-ing"
fdr -i concat/Searchlight_n17_pMap_${maskName}_${modelName} -m ${mDir}/3mm_${maskName} -q $FDRRATE -a corrected/fdrP_${maskName}_${modelName}
#unfortunately, output is p. so to convert to 1-p, do this:
fslmaths corrected/fdrP_${maskName}_${modelName} -mul -1 -add 1 -mas ${mDir}/3mm_${maskName} corrected/fdrP_${maskName}_${modelName}
fslmaths corrected/fdrP_${maskName}_${modelName} -thr .95 -bin corrected/tmp
fslmaths corrected/tmp -Tmean corrected/avg_fdrP_${maskName}_${modelName}
rm corrected/tmp.nii.gz

fdr -i zscored/p_${maskName}_${modelName} -m ${mDir}/3mm_${maskName} -q $FDRRATE -a corrected/fdr_${maskName}_${modelName}
#unfortunately, output is p. so to convert to 1-p, do this:
fslmaths corrected/fdr_${maskName}_${modelName} -mul -1 -add 1 -mas ${mDir}/3mm_${maskName} corrected/fdr_${maskName}_${modelName}
fslmaths corrected/fdr_${maskName}_${modelName} -thr .95 -bin corrected/tmp
fslmaths corrected/tmp -Tmean corrected/avg_fdr_${maskName}_${modelName}
rm corrected/tmp.nii.gz


echo "randomise z values to confirm"
#randomise -i zscored/${maskName}_${modelName} -o corrected/n${n}_${maskName}_${modelName} -1 -m ${mDir}/3mm_${maskName} -x -T --uncorrp -n ${n}

echo "IFG test"
randomise -i zscored/${maskName}_${modelName} -o corrected/n${n}_lIFG_${modelName} -1 -m ${mDir}/left_IFG -x -T --uncorrp -n ${n}

echo "TC test"
randomise -i zscored/${maskName}_${modelName} -o corrected/n${n}_lTC_${modelName} -1 -m ${mDir}/left_TC -x -T --uncorrp -n ${n}
 
:<<same
#this outputs the SAME result as above. no need to do.
echo "randomise z values to confirm"
randomise -i raw/LMVPA_pMap_${maskName}_${modelName} -o corrected/n${n}P_${maskName}_${modelName} -1 -m ${mDir}/3mm_${maskName} -x -T --uncorrp -n ${n}
same
done

done



