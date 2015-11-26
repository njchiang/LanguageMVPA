function [] = ROI_LSA()

% LMVPA_fMRI

% this script performs region of interest analysis on fMRI data.
% Written by Jeff Chiang
% Based on Cai Wingfield 5-2010, 6-2010, 7-2010, 8-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

try
    toolboxRoot = 'Z:\GitHub\LanguageMVPA\multivariate\matlab'; addpath(genpath(toolboxRoot));
catch
    toolboxRoot = '~/GitHub/LanguageMVPA/multivariate/matlab'; addpath(genpath(toolboxRoot));
end
cd Z:\fmri\LanguageMVPA
userOptions = defineUserOptions_LSA;
userOptions.analysisName='ROI_LSA_n16';
gotoDir(userOptions.rootPath);
analysisType='RSA';
userOptions.analysisType=analysisType;
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
% %%%%%%%%%%%%%%%%%%%%%%
% % load data and masks
% % fullBrainVols = fMRIDataPreparation('SPM', userOptions);
fullBrainVols=fMRIDataPreparation(betaCorrespondence_LSA(), userOptions);
binaryMasks_nS = fMRIMaskPreparation(userOptions);
responsePatterns = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence_LSA(), userOptions);
% clear fullBrainVols binaryMasks_nS
load('Z:\fmri\LanguageMVPA\ImageData\ROI_LSA_responsePatterns.mat')

%
% %% RSA %%
% %%%%%%%%%%%%%%%%%%%%%
% %% RDM calculation %%
% %%%%%%%%%%%%%%%%%%%%%
% disp(['Starting RSA analysis']);
sRDMs = constructRDMs(responsePatterns, betaCorrespondence_LSA(), userOptions);

RDMs=averageRDMs_subjectSession(sRDMs, 'session', 'subject');
figureRDMs(RDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1));

LModels=constructModelRDMs(modelRDMs_LSA_L(),userOptions);
PModels=constructModelRDMs(modelRDMs_LSA_P(),userOptions);
% Models=constructModelRDMs(modelRDMs_domain(),userOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LRDMs=RDMs;
PRDMs=RDMs;
for i = 1:length(userOptions.maskNames)
    LRDMs(i).RDM=RDMs(i).RDM(1:32, 1:32);
    PRDMs(i).RDM=RDMs(i).RDM(33:64, 33:64);
end



figureRDMs(LRDMs, userOptions, struct('fileName', 'Lang_RoIRDMs', 'figureNumber', 1));
figureRDMs(PRDMs, userOptions, struct('fileName', 'Pic_RoIRDMs', 'figureNumber', 1));

figureRDMs(RDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1));
figureRDMs(LModels, userOptions, struct('fileName', 'LModelRDMs', 'figureNumber', 2));
figureRDMs(PModels, userOptions, struct('fileName', 'PModelRDMs', 'figureNumber', 2));

%      MDSConditions(sentRDMs, userOptions);
%     dendrogramConditions(sentRDMs, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% relationship amongst multiple RDMs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pairwiseCorrelateRDMs({RDMs, Models}, userOptions);
% pairwiseCorrelateRDMs({sentRDMs, sentModels}, userOptions);
MDSRDMs({RDMs, Models}, userOptions);
MDSRDMs({RDMs, LModels}, userOptions);
MDSRDMs({RDMs, PModels}, userOptions);
%     MDSRDMs({sentRDMs, sentModels}, userOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistical inference %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
userOptions.RDMcorrelationType='Kendall_taua';
% userOptions.RDMcorrelationType='Spearman';
userOptions.RDMrelatednessTest = 'subjectRFXsignedRank';
userOptions.RDMrelatednessThreshold = 0.05;
userOptions.figureIndex = [10 11];
userOptions.RDMrelatednessMultipleTesting = 'FDR';
userOptions.candRDMdifferencesTest = 'subjectRFXsignedRank';
userOptions.candRDMdifferencesThreshold = 0.05;
userOptions.candRDMdifferencesMultipleTesting = 'FDR';
userOptions.significanceTestPermutations = 10000;
userOptions.nResamplings = 10000;
userOptions.nRandomisations=10000;
userOptions.nBootstrap=10000;
userOptions.plotpValues='=';
userOptions.saveFiguresFig = false;
userOptions.saveFiguresPDF=false;
% need to actually separate... i don't really know why
for i = 1: size(sRDMs,1)
    aRDMs{i}=sRDMs(i,:);
end
for i=1:numel(LModels)
    lFmodels{i}=LModels(i);
    pFmodels{i}=PModels(i);
end
% for i=1:numel(Models)
%     models{i}=Models(i);
% end
clear i
close all


for i = 1:numel(aRDMs)
    userOptions.figureIndex = [4*(i-1)+1 (4*(i-1))+2];
    
    roiName=userOptions.maskNames{i};
    disp(['Processing ' roiName])
    userOptions.figure1filename = [userOptions.analysisName '_' roiName '_L_barGraph'];
    userOptions.figure2filename = [userOptions.analysisName '_' roiName '_L_Pvals'];
    %         switched so we match ROIs to model
    %             stats_p_r.(modelName)=compareRefRDM2candRDMs(models, RDMs(roiIndex), userOptions);
    userOptions.resultsPath = fullfile(userOptions.rootPath, 'Statistics/');
    stats_p_r=compareRefRDM2candRDMs(aRDMs{i}, lFmodels, userOptions);
    model2ROIstruct_mfx{i,1}=stats_p_r;
    model2ROIstruct_mfx{i,1}.roiName=roiName;
    userOptions.figureIndex = [4*(i-1)+3 (4*(i-1))+4];
    
    userOptions.figure1filename = [userOptions.analysisName '_' roiName '_P_barGraph'];
    userOptions.figure2filename = [userOptions.analysisName '_' roiName '_P_Pvals'];
    stats_p_r=compareRefRDM2candRDMs(aRDMs{i}, pFmodels, userOptions);
    model2ROIstruct_mfx{i,2}=stats_p_r;
    model2ROIstruct_mfx{i,2}.roiName=roiName;
end

save([userOptions.rootPath '/Statistics/' userOptions.analysisName '_ROI.mat'], '-v7.3')
% save([userOptions.analysisName '_ROI_NaN.mat'], '-v7.3')


end


