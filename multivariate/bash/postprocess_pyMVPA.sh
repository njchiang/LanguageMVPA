#!/bin/bash
:<<doc
registers stat maps to standard space
runs randomise
doc
#standard=${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
projectDir=/Volumes/fmri/LanguageMVPA
#projectDir=/media/sf_fmri/LanguageMVPA
refImage=${projectDir}/fnirt/MNI152_T1_3mm_brain.nii.gz
refMask=${projectDir}/fnirt/MNI152_T1_3mm_brain_mask.nii.gz
targetDir=${projectDir}/Maps/PyMVPA

cd ${targetDir}
#move raw outputs
mkdir raw
mkdir std
mv *.nii.gz raw
cd raw
for model in ActPass_lang ActPass_pic #cross_verb_L2P cross_verb_P2L
do
for indMap in `ls | grep _${model}_cvsl.nii.gz`
do
sub=`echo ${indMap} | cut -d '_' -f1`
echo ${sub}
fslmaths ${indMap} -bin tmp.nii.gz
fslmaths tmp.nii.gz -mul 50 tmp.nii.gz
fslmaths ${indMap} -mul 100 -sub tmp.nii.gz rnd_${indMap}
regMat=${projectDir}/registration/${sub}_example_func2standard.mat

flirt -in rnd_${indMap} -out ../std/std_${indMap} -ref ${refImage} -applyxfm -init ${regMat}
rm rnd_${indMap}
done
fslmerge -t ../Group_${model}.nii.gz ../std/std*_${model}*nii.gz
randomise -i ../Group_${model}.nii.gz -o ../n1000_${model} -1 -T --uncorrp -n 1000 -m ${projectDir}/masks/3mm_grayMatter.hdr
done

