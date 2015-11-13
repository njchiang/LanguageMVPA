function [] = ROI_LSA_across()

% LMVPA_fMRI
%preprocessing was all done, so just go straight to loading RDMs!
% all comparisons will be done on sentence stimuli ONLY
% this script performs region of interest analysis on fMRI data.
% Written by Jeff Chiang
% Based on Cai Wingfield 5-2010, 6-2010, 7-2010, 8-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

toolboxRoot = '~/Box Sync/UCLA/Research/LanguageMVPA/code'; addpath(genpath(toolboxRoot));
% cd /space/raid5/data/monti/Analysis/LanguageMVPA/RSA
cd Z:\fmri\LanguageMVPA
userOptions = defineUserOptions_LSA();
% userOptions = defineUserOptions_Syntax();
gotoDir(userOptions.rootPath);
analysisType='SVM';
userOptions.analysisType=analysisType;
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
% %%%%%%%%%%%%%%%%%%%%%%
% % load data and masks
% % fullBrainVols = fMRIDataPreparation('SPM', userOptions);
% fullBrainVols=fMRIDataPreparation(betaCorrespondence_Semantics(), userOptions);
% binaryMasks_nS = fMRIMaskPreparation(userOptions);
% load(['ImageData/' userOptions.analysisName '_ImageData.mat']);
% load(['ImageData/LSA_masks.mat']);
userOptions.analysisName='SVM_LSA';

% responsePatterns = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence_LSA(), userOptions);
clear fullBrainVols binaryMasks_nS

load(['ImageData/ROI_LSA_responsePatterns.mat']);
% don't need to separate before correlation, they will be the same rank
% order after separating
models=makeLabels_LSA;
% models of interest: 2 3 4 7 8 9
runs=ones(32,1);
nPerm=1000;
svmKernel=2;

% syntax-- within
syntaxLabel=models(12).label(1:32);
syntaxPermutedLabels=permuteLabels(syntaxLabel, runs, nPerm);
actpassLabel=models(13).label(1:32);
actpassLabels=permuteLabels(actpassLabel, runs, nPerm);
relcanLabel=models(14).label(1:32);
relcanPermutedLabels=permuteLabels(relcanLabel, runs, nPerm);
complexLabel=models(15).label(1:32);
complexPermutedLabels=permuteLabels(complexLabel, runs, nPerm);
larpcLabel=models(16).label(1:32);
larpcPermutedLabels=permuteLabels(larpcLabel(larpcLabel>0), runs(larpcLabel>0), nPerm);

verbLabel=models(2).label(1:32);
verbPermutedLabels=permuteLabels(verbLabel, runs, nPerm);
animFamLabel=models(5).label(1:32);
animFamPermutedLabels=permuteLabels(animFamLabel(animFamLabel>0), runs(animFamLabel>0), nPerm);
inAnimFamLabel=models(6).label(1:32);
inAnimFamPermutedLabels=permuteLabels(inAnimFamLabel(inAnimFamLabel>0), runs(inAnimFamLabel>0), nPerm);
famLabel=models(3).label(1:32);
famPermutedLabels=permuteLabels(famLabel, runs, nPerm);
animLabel=models(4).label(1:32);
animPermutedLabels=permuteLabels(animLabel, runs, nPerm);

for m = 1:length(userOptions.maskNames)-1
    currMask=userOptions.maskNames{m};
    for s = 1:length(userOptions.subjectNames)
        currSub=userOptions.subjectNames{s};
        currPat=zscore(responsePatterns.(currMask).(currSub));
        % zscore spatially, THEN transpose
                % separate the correct voxels
[l2l(s,m).Verb, l2p(s,m).Verb, p2p(s,m).Verb, p2l(s,m).Verb, lVerb, pVerb] = run_SVM_LSA_LMVPA(verbLabel, verbPermutedLabels, currPat, runs, svmKernel)   

        sentPats=currPat(:,1:32)';        
        lPats=sentPats(verbLabel>0,:);
        picPats=currPat(:,33:64)';
        pPats=picPats(verbLabel>0,:);
        lVerb(s,m)=optimizeSVM(verbPermutedLabels(:,1), lPats, svmKernel);
        if svmKernel==2
            opts=['-s 0 -t ' num2str(svmKernel) ' -c ' num2str(lVerb(s,m).best_C) ' -g ' num2str(lVerb(s,m).best_gamma)];
        else
            opts=['-s 1 -t ' num2str(svmKernel) ' -n ' num2str(lVerb(s,m).best_C) ];
        end
        lVerbStruct(s,m)=libsvmtrain(verbPermutedLabels(:,1), lPats, opts);
        [l2l(s,m).Verb.acc, l2l(s,m).Verb.dist, l2l(s,m).Verb.thresh]=CVpermutationSVM(lPats, verbPermutedLabels, runs, opts);
        [l2p(s,m).Verb.acc, l2p(s,m).Verb.dist, l2p(s,m).Verb.thresh] = permutationSVM(lVerbStruct, pPats, verbPermutedLabels, runs);
  
        pVerb(s,m)=optimizeSVM(verbPermutedLabels(;,1), pPats, svmKernel);
        if svmKernel==2
            opts=['-s 0 -t ' num2str(svmKernel) ' -c ' num2str(pVerb(s,m).best_C) ' -g ' num2str(pVerb(s,m).best_gamma)];
        else
            opts=['-s 1 -t ' num2str(svmKernel) ' -n ' num2str(pVerb(s,m).best_C) ];
        end
        pVerbStruct(s,m)=libsvmtrain(verbPermutedLabels(:,1), pPats, opts);
        [p2p(s,m).Verb.acc, p2p(s,m).Verb.dist, p2p(s,m).Verb.thresh]=CVpermutationSVM(pPats, verbPermutedLabels, runs, opts);
        [p2l(s,m).Verb.acc, p2l(s,m).Verb.dist, p2l(s,m).Verb.thresh] = permutationSVM(pVerbStruct, sentPats, verbPermutedLabels, runs);

    end
end

save('Statistics/ROI_SVM.mat', 'l2l', 'l2p', 'p2p', 'p2l', 'lAnimStruct', 'pAnimStruct', 'lVerbStruct', 'pVerbStruct', 'lFamStruct', 'pFamStruct' )

