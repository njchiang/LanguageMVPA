#!/bin/bash
:<<doc
This script gets multi-number measures per subject, concatenates them properly
and runs randomise.
doc

model=${3}
mask=${2}
chance=${4}
t=${5}
projectDir=/Volumes/fmri/LanguageMVPA
#projectDir=/Volumes/JEFF/UCLA/LMVPA
desDir=/Users/njchiang/GitHub/LanguageMVPA/multivariate/bash
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

fslmerge -t ../${mask}_${model}_Group.nii.gz \
../std/std_LMVPA001_${mask}_${model}.nii.gz \
../std/std_LMVPA002_${mask}_${model}.nii.gz \
../std/std_LMVPA003_${mask}_${model}.nii.gz \
../std/std_LMVPA005_${mask}_${model}.nii.gz \
../std/std_LMVPA006_${mask}_${model}.nii.gz \
../std/std_LMVPA007_${mask}_${model}.nii.gz \
../std/std_LMVPA008_${mask}_${model}.nii.gz \
../std/std_LMVPA009_${mask}_${model}.nii.gz \
../std/std_LMVPA010_${mask}_${model}.nii.gz \
../std/std_LMVPA011_${mask}_${model}.nii.gz \
../std/std_LMVPA013_${mask}_${model}.nii.gz \
../std/std_LMVPA014_${mask}_${model}.nii.gz \
../std/std_LMVPA015_${mask}_${model}.nii.gz \
../std/std_LMVPA016_${mask}_${model}.nii.gz \
../std/std_LMVPA017_${mask}_${model}.nii.gz \
../std/std_LMVPA018_${mask}_${model}.nii.gz \
../std/std_LMVPA019_${mask}_${model}.nii.gz 

randomise -i ../${mask}_${model}_Group.nii.gz -o ../n1000_${mask}_${model} \
	-v 5 -d $desDir/${t}mfx_full_design.mat -t $desDir/${t}mfx_full_design.con \
	-T --uncorrp -n 1000 -m ${projectDir}/masks/3mm_grayMatter
