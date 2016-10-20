#!/bin/bash
:<<doc
registers each run to the example_func!
doc

projectDir=/Volumes/fmri/LanguageMVPA/data

# this part does regressed data
# 001 002 003 005 006 007 008 009 010 011 013 014 015 016 017 018 019
:<<func

for i in 001 002 003 005 006 007 008 009 010 011 013 014 015 016 017 018 019
do
	s=LMVPA${i}
	echo ${s}
	cd ${projectDir}/${s}/func/extra
	for scan in ${s}_*_mc.nii.gz
	do
		r=`echo ${scan} | cut -d '_' -f2`
		scp njchiang@psych.ucla.edu:/space/raid5/data/monti/Analysis/LanguageMVPA/${s}/${s}_${r}_mc.feat/stats/res4d.nii.gz ${s}_${r}_mc_reg.nii.gz
		if [[ "$r" == "Run1" ]] 
		then
			cp ${scan} ../${s}_${r}.nii.gz
		else
			flirt -in ${scan} -ref example_func.nii.gz -out ../${s}_${r}.nii.gz -applyxfm -init ../../reg/${r}_to_1.mat
		fi
	done
done	
for sub in LMVPA001 LMVPA002 LMVPA003 LMVPA005 LMVPA006 LMVPA007 LMVPA008 LMVPA009 LMVPA010 LMVPA011 LMVPA013 LMVPA014 LMVPA015 LMVPA016 LMVPA017 LMVPA018 LMVPA019 LMVPA020
do
  cd ${projectDir}/${sub}
  for scan in ${sub}_Run?.nii.gz
  do
    echo ${scan}
    fName=`echo ${scan} | cut -d '.' -f1`
    run=`echo ${fName} | cut -d '_' -f2`
    mkdir ${fName}_mc
    fslroi ${scan} ${fName}_mc 4 -1
    mcflirt -in ${fName}_mc -out ${fName}_mc -mats -plots -reffile ${run}_preproc.feat/example_func -rmsrel -rmsabs -spline_final
    mv -f ${fName}_mc.mat ${fName}_mc.par ${fName}_mc_abs.rms ${fName}_mc_abs_mean.rms ${fName}_mc_rel.rms ${fName}_mc_rel_mean.rms ${fName}_mc

  done
done
func

# The new and improved preprocessing for an encoding model analysis. Minimal preprocessing is done. Really just deleting volumes, motion correction and slicetiming.

projectDir=~/data/fmri/LanguageMVPA/data
for i in 001 002 003 005 006 007 008 009 010 011 013 014 015 016 017 018 019
do
	s=LMVPA${i}
	echo ${s}
	cd ${projectDir}/${s}/func
	rm extra/*_preprocessed.nii.gz
	for scan in ${s}_Run?.tsv
	do
		r=`echo ${scan} | cut -d '.' -f1 | cut -d '_' -f2`
		echo ${s} ${r}
		echo "cleaning..."
		rm ${s}_${r}.nii.gz # remove whatever version i have now
		mv ${s}_${r}_reg.nii.gz ${s}_${r}_resid.nii.gz
		mv extra/${s}_${r}_mc_reg.nii.gz extra/${s}_${r}_mc_resid.nii.gz
		echo "slicetiming..."
		slicetimer -i extra/${s}_${r}_mc -o extra/${s}_${r}_mc_st -r 3 --ocustom=${projectDir}/SliceOrdering.txt
		echo "registering"
		if [[ "$r" == "Run1" ]] 
		then
			cp extra/${s}_${r}_mc_st.nii.gz ${s}_${r}.nii.gz
		else
			flirt -in extra/${s}_${r}_mc_st.nii.gz -ref extra/example_func.nii.gz -out ${s}_${r}.nii.gz -applyxfm -init ../reg/${r}_to_1.mat
		fi
		
	done
done	