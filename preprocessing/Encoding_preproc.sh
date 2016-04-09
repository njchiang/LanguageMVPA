#!/bin/bash
:<<doc
The new and improved preprocessing for an encoding model analysis. Minimal preprocessing is done. Really just deleting volumes and preprocessing.
doc

projectDir=/space/raid5/data/monti/Analysis/LanguageMVPA

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
