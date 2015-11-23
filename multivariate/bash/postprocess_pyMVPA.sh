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
:<<com
for model in syntax_lang syntax_pic anim_verb_lang anim_verb_pic \
	inanim_verb_lang inanim_verb_pic 
do
	for indMap in `ls | grep grayMatter_${model}_concat_cvsl.nii.gz`
	do
		sub=`echo ${indMap} | cut -d '_' -f1`
		meanMap=${sub}_${model}
		fslmaths ${indMap} -Tmean ${meanMap}
		echo ${sub}
		fslmaths ${projectDir}/masks/sub/${sub}_grayMatter.nii.gz -bin tmp.nii.gz
		fslmaths tmp.nii.gz -mul 25 tmp.nii.gz
		fslmaths ${meanMap} -mul 100 -sub tmp.nii.gz rnd_${meanMap}
		regMat=${projectDir}/registration/${sub}_example_func2standard.mat

		flirt -in rnd_${meanMap} -out ../std/std_${meanMap} -ref ${refImage} \
			-applyxfm -init ${regMat}
		rm rnd_${meanMap}.nii.gz
	done
fslmerge -t ../Group_${model}.nii.gz ../std/std*_${model}*nii.gz
randomise -i ../Group_${model}.nii.gz -o ../n1000_${model} -1 -T --uncorrp -n 1000 -m ${projectDir}/masks/3mm_grayMatter.hdr
done

for model in ActPass_lang ActPass_pic RelCan_lang RelCan_pic anim_ActPass_lang \
	anim_ActPass_pic inanim_ActPass_lang inanim_ActPass_pic anim_RelCan_lang \
	anim_RelCan_pic inanim_RelCan_lang inanim_RelCan_pic syntax_pic
do
	for indMap in `ls | grep grayMatter_${model}_concat_cvsl.nii.gz`
	do
		sub=`echo ${indMap} | cut -d '_' -f1`
		meanMap=${sub}_${model}
		fslmaths ${indMap} -Tmean ${meanMap}
		echo ${sub}
		fslmaths ${projectDir}/masks/sub/${sub}_grayMatter.nii.gz -bin tmp.nii.gz
		fslmaths tmp.nii.gz -mul 50 tmp.nii.gz
		fslmaths ${meanMap} -mul 100 -sub tmp.nii.gz rnd_${meanMap}
		regMat=${projectDir}/registration/${sub}_example_func2standard.mat

		flirt -in rnd_${meanMap} -out ../std/std_${meanMap} -ref ${refImage} \
			-applyxfm -init ${regMat}
		rm rnd_${meanMap}.nii.gz
	done
fslmerge -t ../Group_${model}.nii.gz ../std/std*_${model}*nii.gz
randomise -i ../Group_${model}.nii.gz -o ../n1000_${model} -1 -T --uncorrp -n 1000 -m ${projectDir}/masks/3mm_grayMatter.hdr
done
com
for model in verb_lang verb_pic
do
	for indMap in `ls | grep grayMatter_${model}_concat_cvsl.nii.gz`
	do
		sub=`echo ${indMap} | cut -d '_' -f1`
		meanMap=${sub}_${model}
		fslmaths ${indMap} -Tmean ${meanMap}
		echo ${sub}
		fslmaths ${projectDir}/masks/sub/${sub}_grayMatter.nii.gz -bin tmp.nii.gz
		fslmaths tmp.nii.gz -mul 12.5 tmp.nii.gz
		fslmaths ${meanMap} -mul 100 -sub tmp.nii.gz rnd_${meanMap}
		regMat=${projectDir}/registration/${sub}_example_func2standard.mat

		flirt -in rnd_${meanMap} -out ../std/std_${meanMap} -ref ${refImage} \
			-applyxfm -init ${regMat}
		rm rnd_${meanMap}.nii.gz
	done
fslmerge -t ../Group_${model}.nii.gz ../std/std*_${model}*nii.gz
randomise -i ../Group_${model}.nii.gz -o ../n1000_${model} -1 -T --uncorrp -n 1000 -m ${projectDir}/masks/3mm_grayMatter.hdr
done

