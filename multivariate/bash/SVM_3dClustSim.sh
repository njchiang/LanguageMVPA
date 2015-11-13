#!/bin/bash
:<<doc
gets ready for and runs randomise

doc

cd /Volumes/fmri/LanguageMVPA/Maps/Syntax/SVM
for d in Group_*nii.gz
	do
		model=`echo ${d} | cut -d '_' -f2 | cut -d '.' -f1`
		echo ${model}
		sh /Volumes/fmri/LanguageMVPA/run_3dClustSim.sh ${model} ${d}
	done

cd /Volumes/fmri/LanguageMVPA/Maps/Semantics/SVM
for d in Group_*nii.gz
	do
		model=`echo ${d} | cut -d '_' -f2 | cut -d '.' -f1`
		echo ${model}
		sh /Volumes/fmri/LanguageMVPA/run_3dClustSim.sh ${model} ${d}
	done

