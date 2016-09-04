#!/bin/bash
:<<doc
newDirName
inputImage
inputMask
doc

newDirName=${1}
mkdir ${newDirName}
echo "Copying..."
cp ${2} ${newDirName}/rawImage.nii.gz
#mask=/media/sf_fmri/LanguageMVPA/masks/3mm_grayMatter.nii.gz
mask=/Volumes/fmri/LanguageMVPA/masks/3mm_grayMatter.nii.gz

cd ${newDirName}

echo "splitting..."
fslsplit rawImage.nii.gz splitImage -t

echo "t-testing..."
3dttest++ -prefix 1Sample.nii.gz -setA splitImage*.nii.gz -toz

echo "estimating smoothness..."
echo `3dFWHMx -mask ${mask} 1Sample.nii.gz` > fwhm.txt

fwhm=`cat fwhm.txt`

echo "3dClustSim..."
x=`fslval rawImage dim1`
y=`fslval rawImage dim2`
z=`fslval rawImage dim3`
px=`fslval rawImage pixdim1`
py=`fslval rawImage pixdim2`
pz=`fslval rawImage pixdim3`

3dClustSim -mask ${mask} -fwhmxyz ${fwhm} -nxyz ${x} ${y} ${z} -dxyz ${px} ${py} ${pz} -nodec -niter 1000 -prefix 1Sample.nii.gz -NN 123

#find corresponding entry... not sure which files to look at.
#1.65 2.33
tCutoff=2.33 #p=.01 uncorrected cutoff
if [ -f 1Sample.nii.gz.NN2_1sided.1D ]
then
tmp=`sed '9q;d' 1Sample.nii.gz.NN2_1sided.1D`
elif [ -f 1Sample.nii.gz.NN2.1D ]
then
tmp=`sed '9q;d' 1Sample.nii.gz.NN2.1D`
else
tmp=`sed '9q;d' 1Sample.nii.gz.NN1.1D`
fi
clustSize=`echo ${tmp} | cut -d ' ' -f3` #p=.05 cluster size

echo "output"
${FSLDIR}/bin/cluster -z 1Sample.nii.gz -t ${tCutoff} --minextent=${clustSize} --othresh=../${newDirName}_p01_c05_image.nii.gz > ../${newDirName}_p01_c05_clusters.txt

echo "cleaning up"
cd ..
rm -r ${newDirName}