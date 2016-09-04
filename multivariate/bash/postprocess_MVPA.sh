#!/bin/bash
:<<doc
registers stat maps to standard space
runs randomise
usage:
sh postprocess_LMVPA.sh PyMVPA grayMatter cross_anim_L2P_ccsl 50
doc

#standard=${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
model=${3}
mask=${2}
chance=${4}
projectDir=/Volumes/fmri/LanguageMVPA
#projectDir=/Volumes/JEFF/UCLA/LMVPA
#projectDir=/media/sf_fmri/LanguageMVPA
refImage=${projectDir}/MNI152_T1_3mm_brain.nii.gz
refMask=${projectDir}/fnirt/MNI152_T1_3mm_brain_mask.nii.gz
targetDir=${projectDir}/Maps/${1}

cd ${targetDir}
#move raw outputs
mkdir raw
mkdir std
mv *_${mask}_${model}.nii.gz raw
cd raw
for indMap in `ls | grep _${mask}_${model}.nii.gz`
do
sub=`echo ${indMap} | cut -d '_' -f1`
echo ${sub}
fslmaths ${projectDir}/data/$sub/masks/${sub}_grayMatter.nii.gz -bin tmp.nii.gz

if [ "$chance" = "0" ]
then
	cp ${indMap} rnd_${indMap}
else
	fslmaths tmp.nii.gz -mul ${chance} tmp.nii.gz
	fslmaths ${indMap} -mul 100 -sub tmp.nii.gz rnd_${indMap}
fi


regMat=${projectDir}/data/$sub/reg/${sub}_example_func2standard.mat

flirt -in rnd_${indMap} -out ../std/std_${indMap} -ref ${refImage} \
	-applyxfm -init ${regMat}
rm rnd_${indMap} tmp.nii.gz
done

fslmerge -t ../${mask}_${model}_Group.nii.gz ../std/std*_${mask}_${model}*nii.gz
randomise -i ../${mask}_${model}_Group.nii.gz -o ../n1000_${mask}_${model} \
	-v 5 -1 -T -x --uncorrp -n 1000 -m ${projectDir}/data/standard/3mm_grayMatter


