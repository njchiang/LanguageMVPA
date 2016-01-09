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
targetDir=${projectDir}/Maps/LSA

cd ${targetDir}
#move raw outputs
mkdir raw
mkdir std
mv *.hdr raw/
mv *.img raw/
cd raw
echo Converting...
for f in *.hdr
do
	fslchfiletype NIFTI_GZ ${f}
done
echo Working...
for h in L P C
do
for m2 in WholeTopics WholeLSA VerbLSA SubjectLSA ObjectLSA NPLSA VerbDetector SyntaxDetector SyntaxComplex AnimDetector
do
	model=${h}${m2}
	for indMap in `ls | grep grayMatter_${model}_rMap.nii.gz`
	do
		sub=`echo ${indMap} | cut -d '_' -f1`
		echo ${sub}
		regMat=${projectDir}/registration/${sub}_example_func2standard.mat

		flirt -in ${indMap} -out ../std/std_${indMap} -ref ${refImage} \
			-applyxfm -init ${regMat}
	done
fslmerge -t ../Group_${model}.nii.gz ../std/std*_${model}*nii.gz
randomise -i ../Group_${model}.nii.gz -o ../n1000_${model} -1 -T --uncorrp -n 1000 -m ${projectDir}/masks/3mm_grayMatter.hdr
done
done

for model in WholeTopics WholeLSA SubjectLSA ObjectLSA NPLSA VerbDetector SyntaxDetector SyntaxComplex AnimDetector
do
	for indMap in `ls | grep grayMatter_${model}_rMap.nii.gz`
	do
		sub=`echo ${indMap} | cut -d '_' -f1`
		echo ${sub}
		regMat=${projectDir}/registration/${sub}_example_func2standard.mat
		flirt -in ${indMap} -out ../std/std_${indMap} -ref ${refImage} \
			-applyxfm -init ${regMat}
	done
fslmerge -t ../Group_${model}.nii.gz ../std/std*_${model}*nii.gz
randomise -i ../Group_${model}.nii.gz -o ../n1000_${model} -1 -T --uncorrp -n 1000 -m ${projectDir}/masks/3mm_grayMatter.hdr
done

