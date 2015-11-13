#!/bin/bash

analysisType=$1
mapsDir=/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/Maps/${analysisType}/raw
outDir=/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/Maps/${analysisType}
masksDir=/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/data/masks
mkdir -p ${outDir}/concat
mkdir -p ${outDir}/nMaps-9mm
cd ${mapsDir}

for mask in `ls ${masksDir}/3mm_grayMatter.hdr`
do

rheaderName=`echo ${mask} | cut -d'.' -f1 | cut -d'/' -f11 | sed "s/3mm_/Searchlight_n17_rMap_/g"`
pheaderName=`echo ${mask} | cut -d'.' -f1 | cut -d'/' -f11 | sed "s/3mm_/Searchlight_n17_pMap_/g"`
nheaderName=`echo ${mask} | cut -d'.' -f1 | cut -d'/' -f11 | sed "s/3mm_/Searchlight_n17_nMap_/g"`

echo ${nheaderName}
cp ${mapsDir}/${nheaderName}_LMVPA001.nii.gz /space/raid5/data/monti/Analysis/LanguageMVPA/RSA/Maps/${analysisType}/nMaps-9mm/${nheaderName}.nii.gz 
for modelName in Language Pictures StimType LvP Syntax AllSemantics Semantics LSemantics PSemantics LvPSemantics
do
echo Merging
echo ${rheaderName}_${modelName}
fslmerge -t ${outDir}/concat/${rheaderName}_${modelName}_concat \
${mapsDir}/${rheaderName}_${modelName}_LMVPA001 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA002 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA003 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA005 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA006 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA007 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA008 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA009 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA010 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA011 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA013 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA014 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA015 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA016 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA017 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA018 \
${mapsDir}/${rheaderName}_${modelName}_LMVPA019

#${mapsDir}/${rheaderName}_${modelName}_LMVPA007 \

fslmaths ${outDir}/concat/${rheaderName}_${modelName}_concat -nan ${outDir}/concat/${rheaderName}_${modelName}
echo ${pheaderName}_${modelName}
fslmerge -t ${outDir}/concat/${pheaderName}_${modelName}_concat \
${mapsDir}/${pheaderName}_${modelName}_LMVPA001 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA002 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA003 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA005 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA006 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA007 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA008 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA009 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA010 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA011 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA013 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA014 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA015 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA016 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA017 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA018 \
${mapsDir}/${pheaderName}_${modelName}_LMVPA019

#${mapsDir}/${pheaderName}_${modelName}_LMVPA007 \
fslmaths ${outDir}/concat/${rheaderName}_${modelName}_concat -nan ${outDir}/concat/${pheaderName}_${modelName}
done
done
