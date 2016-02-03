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
echo Working...
for m2 in WholeTopics WholeLSA VerbLSA SubjectLSA ObjectLSA NPLSA VerbDetector SyntaxDetector SyntaxComplex AnimDetector
do
	echo ${m2}
	sh ~/GitHub/LanguageMVPA/multivariate/bash/postprocess_MVPA.sh LSA grayMatter ${m2}_rMap 0
	
for h in L P C
do 
	model=${h}${m2}_rMap
	echo ${model}
	sh ~/GitHub/LanguageMVPA/multivariate/bash/postprocess_MVPA.sh LSA grayMatter ${model} 0
done
done


