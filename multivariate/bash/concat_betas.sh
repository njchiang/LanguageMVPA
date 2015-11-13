#!/bin/bash
:<<doc

doc

cd /Volumes/fmri/LanguageMVPA/betas
for sub in LMVPA001 LMVPA002 LMVPA003 LMVPA005 LMVPA006 LMVPA007 LMVPA008 LMVPA009 LMVPA010 LMVPA011 LMVPA013 LMVPA014 LMVPA015 LMVPA016 LMVPA017 LMVPA018 LMVPA019
do
	echo ${sub}

	
fslmerge -t ${sub}_LSA_Series.nii.gz \
${sub}_s_Touch_Act_Can.nii.gz ${sub}_s_Touch_Act_Rel.nii.gz \
${sub}_s_Touch_Pass_Can.nii.gz ${sub}_s_Touch_Pass_Rel.nii.gz \
${sub}_s_Stretch_Act_Can.nii.gz ${sub}_s_Stretch_Act_Rel.nii.gz \
${sub}_s_Stretch_Pass_Can.nii.gz ${sub}_s_Stretch_Pass_Rel.nii.gz \
${sub}_s_Light_Act_Can.nii.gz ${sub}_s_Light_Act_Rel.nii.gz \
${sub}_s_Light_Pass_Can.nii.gz ${sub}_s_Light_Pass_Rel.nii.gz \
${sub}_s_Kiss_Act_Can.nii.gz ${sub}_s_Kiss_Act_Rel.nii.gz \
${sub}_s_Kiss_Pass_Can.nii.gz ${sub}_s_Kiss_Pass_Rel.nii.gz \
${sub}_s_Kick_Act_Can.nii.gz ${sub}_s_Kick_Act_Rel.nii.gz \
${sub}_s_Kick_Pass_Can.nii.gz ${sub}_s_Kick_Pass_Rel.nii.gz \
${sub}_s_Hit_Act_Can.nii.gz ${sub}_s_Hit_Act_Rel.nii.gz \
${sub}_s_Hit_Pass_Can.nii.gz ${sub}_s_Hit_Pass_Rel.nii.gz \
${sub}_s_Crush_Act_Can.nii.gz ${sub}_s_Crush_Act_Rel.nii.gz \
${sub}_s_Crush_Pass_Can.nii.gz ${sub}_s_Crush_Pass_Rel.nii.gz \
${sub}_s_Console_Act_Can.nii.gz ${sub}_s_Console_Act_Rel.nii.gz \
${sub}_s_Console_Pass_Can.nii.gz ${sub}_s_Console_Pass_Rel.nii.gz \
${sub}_l_Touch_Act_Can.nii.gz ${sub}_l_Touch_Act_Rel.nii.gz \
${sub}_l_Touch_Pass_Can.nii.gz ${sub}_l_Touch_Pass_Rel.nii.gz \
${sub}_l_Stretch_Act_Can.nii.gz ${sub}_l_Stretch_Act_Rel.nii.gz \
${sub}_l_Stretch_Pass_Can.nii.gz ${sub}_l_Stretch_Pass_Rel.nii.gz \
${sub}_l_Light_Act_Can.nii.gz ${sub}_l_Light_Act_Rel.nii.gz \
${sub}_l_Light_Pass_Can.nii.gz ${sub}_l_Light_Pass_Rel.nii.gz \
${sub}_l_Kiss_Act_Can.nii.gz ${sub}_l_Kiss_Act_Rel.nii.gz \
${sub}_l_Kiss_Pass_Can.nii.gz ${sub}_l_Kiss_Pass_Rel.nii.gz \
${sub}_l_Kick_Act_Can.nii.gz ${sub}_l_Kick_Act_Rel.nii.gz \
${sub}_l_Kick_Pass_Can.nii.gz ${sub}_l_Kick_Pass_Rel.nii.gz \
${sub}_l_Hit_Act_Can.nii.gz ${sub}_l_Hit_Act_Rel.nii.gz \
${sub}_l_Hit_Pass_Can.nii.gz ${sub}_l_Hit_Pass_Rel.nii.gz \
${sub}_l_Crush_Act_Can.nii.gz ${sub}_l_Crush_Act_Rel.nii.gz \
${sub}_l_Crush_Pass_Can.nii.gz ${sub}_l_Crush_Pass_Rel.nii.gz \
${sub}_l_Console_Act_Can.nii.gz ${sub}_l_Console_Act_Rel.nii.gz \
${sub}_l_Console_Pass_Can.nii.gz ${sub}_l_Console_Pass_Rel.nii.gz
done

