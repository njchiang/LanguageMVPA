function [] = ROI_LSA()

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
userOptions = defineUserOptions_Group();
% userOptions = defineUserOptions_Syntax();
% userOptions.analysisName='ROI_LSA';
gotoDir(userOptions.rootPath);
analysisType='RSA';
userOptions.analysisType=analysisType;
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
% %%%%%%%%%%%%%%%%%%%%%%
% % load data and masks
% % fullBrainVols = fMRIDataPreparation('SPM', userOptions);
% fullBrainVols=fMRIDataPreparation(betaCorrespondence_Semantics(), userOptions);
fullBrainVols=fMRIDataPreparation(betaCorrespondence_Group(), userOptions);
binaryMasks_nS = fMRIMaskPreparation(userOptions);
% load(['ImageData/' userOptions.analysisName '_ImageData.mat']);
% load(['ImageData/ROI_masks.mat']);
responsePatterns = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence_Group(), userOptions);
% clear fullBrainVols binaryMasks_nS

% load(['ImageData/' userOptions.analysisName '_responsePatterns.mat']);
% don't need to separate before correlation, they will be the same rank
% order after separating
% zResponsePatternsSpace=responsePatterns;
% zResponsePatternsTime=responsePatterns;

semanticsResponsePatterns=responsePatterns;
syntaxResponsePatterns=responsePatterns;
animResponsePatterns=responsePatterns;
for m = 1:length(userOptions.maskNames)
    currMask=userOptions.maskNames{m};
    for s = 1:length(userOptions.subjectNames)
        currSub=userOptions.subjectNames{s};
        currPat=responsePatterns.(currMask).(currSub);
        semanticsResponsePatterns.(currMask).(currSub)=currPat(:,1:16);
        sSemanticsResponsePatterns.(currMask).(currSub)=currPat(:,1:8);
        pSemanticsResponsePatterns.(currMask).(currSub)=currPat(:,9:16); 
        syntaxResponsePatterns.(currMask).(currSub)=currPat(:,17:24);
        sSyntaxResponsePatterns.(currMask).(currSub)=currPat(:,17:20);
        pSyntaxResponsePatterns.(currMask).(currSub)=currPat(:,21:24);
        animResponsePatterns.(currMask).(currSub)=currPat(:,25:28);
        sAnimResponsePatterns.(currMask).(currSub)=currPat(:,25:26);
        pInanimResponsePatterns.(currMask).(currSub)=currPat(:,27:28);
    end
end

% 
% %% RSA %%
% %%%%%%%%%%%%%%%%%%%%%
% %% RDM calculation %%
% %%%%%%%%%%%%%%%%%%%%%
% disp(['Starting RSA analysis']);
semanticsRDMs = constructRDMs(semanticsResponsePatterns, betaCorrespondence_Group(), userOptions);
semanticsRDMs = averageRDMs_subjectSession(semanticsRDMs, 'session', 'subject');
sSemanticsRDMs = constructRDMs(sSemanticsResponsePatterns, betaCorrespondence_Group(), userOptions);
sSemanticsRDMs = averageRDMs_subjectSession(sSemanticsRDMs, 'session', 'subject');
pSemanticsRDMs = constructRDMs(pSemanticsResponsePatterns, betaCorrespondence_Group(), userOptions);
pSemanticsRDMs = averageRDMs_subjectSession(pSemanticsRDMs, 'session', 'subject');

syntaxRDMs = constructRDMs(syntaxResponsePatterns, betaCorrespondence_Group(), userOptions);
syntaxRDMs = averageRDMs_subjectSession(syntaxRDMs, 'session', 'subject');
sSyntaxRDMs = constructRDMs(sSyntaxResponsePatterns, betaCorrespondence_Group(), userOptions);
sSyntaxRDMs = averageRDMs_subjectSession(sSyntaxRDMs, 'session', 'subject');
pSyntaxRDMs = constructRDMs(pSyntaxResponsePatterns, betaCorrespondence_Group(), userOptions);
pSyntaxRDMs = averageRDMs_subjectSession(pSyntaxRDMs, 'session', 'subject');
animRDMs = constructRDMs(animResponsePatterns, betaCorrespondence_Group(), userOptions);
animRDMs = averageRDMs_subjectSession(animRDMs, 'session', 'subject');
sAnimRDMs = constructRDMs(sAnimResponsePatterns, betaCorrespondence_Group(), userOptions);
sAnimRDMs = averageRDMs_subjectSession(sAnimRDMs, 'session', 'subject');
pAnimRDMs = constructRDMs(pInanimResponsePatterns, betaCorrespondence_Group(), userOptions);
pAnimRDMs = averageRDMs_subjectSession(pAnimRDMs, 'session', 'subject');

figureRDMs(semanticsRDMs, userOptions, struct('fileName', 'Semantics_RoIRDMs', 'figureNumber', 1));
figureRDMs(sSemanticsRDMs, userOptions, struct('fileName', 'SSemantics_RoIRDMs', 'figureNumber', 2));
figureRDMs(pSemanticsRDMs, userOptions, struct('fileName', 'PSemantics_RoIRDMs', 'figureNumber', 3));
figureRDMs(syntaxRDMs, userOptions, struct('fileName', 'Syntax_RoIRDMs', 'figureNumber', 4));
figureRDMs(sSyntaxRDMs, userOptions, struct('fileName', 'SSyntax_RoIRDMs', 'figureNumber', 5));
figureRDMs(pSyntaxRDMs, userOptions, struct('fileName', 'PSyntax_RoIRDMs', 'figureNumber', 6));
figureRDMs(animRDMs, userOptions, struct('fileName', 'Anim_RoIRDMs', 'figureNumber', 7));
figureRDMs(sAnimRDMs, userOptions, struct('fileName', 'SAnim_RoIRDMs', 'figureNumber', 8));
figureRDMs(pAnimRDMs, userOptions, struct('fileName', 'PAnim_RoIRDMs', 'figureNumber', 9));


% RDMs = constructRDMs(responsePatterns, betaCorrespondence_Semantics(), userOptions);


load(['RDMs/' userOptions.analysisName '_RDMs.mat']);

% sRDMs = averageRDMs_subjectSession(RDMs, 'session');
% RDMs = averageRDMs_subjectSession(RDMs, 'session', 'subject');

Models=constructModelRDMs(modelRDMs_ROI_layer(),userOptions);
% Models=constructModelRDMs(modelRDMs_fullSemantics(),userOptions);
% LModels=constructModelRDMs(modelRDMs_ROI_LSyntax(),userOptions);
% PModels=constructModelRDMs(modelRDMs_ROI_PSyntax(),userOptions);
LModels=constructModelRDMs(modelRDMs_LSA_L(),userOptions);
PModels=constructModelRDMs(modelRDMs_LSA_P(),userOptions);
% Models=constructModelRDMs(modelRDMs_domain(),userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



     MDSConditions(sSemanticsRDMs, userOptions);
     MDSConditions(sSyntaxRDMs, userOptions);
allConds=userOptions.conditionLabels;
allCols=userOptions.conditionColours;
userOptions.conditionLabels=allConds{1:8};
userOptions.conditionColors=allCols(1:8,:);
    dendrogramConditions(sSemanticsRDMs, userOptions);

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
%     lRDMs{i}=sRDMs(i,:);
%     pRDMs{i}=sRDMs(i,:);
%     for j = 1:length(userOptions.subjectNames)
%         lRDMs{i}(j).RDM=sRDMs(i,j).RDM(1:32,1:32);
%         pRDMs{i}(j).RDM=sRDMs(i,j).RDM(33:64, 33:64);
%     end
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
    userOptions.figureIndex = [2*(i-1)+1 (2*(i-1))+2];

    roiName=userOptions.maskNames{i};
    disp(['Processing ' roiName])
%     %set other mode model RDM to empirical
    userOptions.figure1filename = [userOptions.analysisName '_' roiName '_barGraph'];
    userOptions.figure2filename = [userOptions.analysisName '_' roiName '_Pvals'];
%     %         switched so we match ROIs to model
%     %             stats_p_r.(modelName)=compareRefRDM2candRDMs(models, RDMs(roiIndex), userOptions);
%     userOptions.resultsPath = fullfile(userOptions.rootPath, 'Statistics/');
    stats_p_r=compareRefRDM2candRDMs(aRDMs{i}, models, userOptions);
    model2ROIstruct_mfx{i,1}=stats_p_r;
    model2ROIstruct_mfx{i,1}.roiName=roiName;

end


%%% separating L and P
 % using NaNs
for i = 1:numel(aRDMs)
    userOptions.figureIndex = [4*(i-1)+1 (4*(i-1))+2];

    roiName=userOptions.maskNames{i};
    disp(['Processing ' roiName])
%     %set other mode model RDM to empirical
%     userOptions.figure1filename = [userOptions.analysisName '_' roiName '_barGraph'];
%     userOptions.figure2filename = [userOptions.analysisName '_' roiName '_Pvals'];
%     %         switched so we match ROIs to model
%     %             stats_p_r.(modelName)=compareRefRDM2candRDMs(models, RDMs(roiIndex), userOptions);
%     userOptions.resultsPath = fullfile(userOptions.rootPath, 'Statistics/');
%     stats_p_r=compareRefRDM2candRDMs(aRDMs{i}, models, userOptions);
%     model2ROIstruct_mfx{i,1}=stats_p_r;
%     model2ROIstruct_mfx{i,1}.roiName=roiName;
%     lFmodels{end}.RDM(1:32,1:32)=RDMs(i).RDM(33:64, 33:64);
%     pFmodels{end}.RDM(33:64, 33:64)=RDMs(i).RDM(1:32, 1:32);
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
% 
% % ActuL sep
% for i=1:numel(LModels)
%     lSmodels{i}=LModels(i);
%     lSmodels{i}.RDM=LModels(i).RDM(1:32,1:32);
%     pSmodels{i}=PModels(i);
%     pSmodels{i}.RDM=PModels(i).RDM(33:64,33:64);
% end
% close all
% 
% for i = 1:numel(aRDMs)
%     userOptions.figureIndex = [4*(i-1)+1 (4*(i-1))+2];
%     maskName=userOptions.maskNames{i};
%     disp(['Processing ' roiName])
%     %set other mode model RDM to empirical
%     lSmodels{end}.RDM=RDMs(i).RDM(33:64, 33:64);
%     pSmodels{end}.RDM=RDMs(i).RDM(1:32, 1:32);
%     userOptions.figure1filename = [userOptions.analysisName '_' roiName '_L_barGraph'];
%     userOptions.figure2filename = [userOptions.analysisName '_' roiName '_L_Pvals'];
%     %         switched so we match ROIs to model
%     %             stats_p_r.(modelName)=compareRefRDM2candRDMs(models, RDMs(roiIndex), userOptions);
%     userOptions.resultsPath = fullfile(userOptions.rootPath, 'Statistics/');
%     stats_p_r=compareRefRDM2candRDMs(lRDMs{i}, lSmodels, userOptions);
%     close all;
%     model2ROIstruct_mfx_sep{i,1}=stats_p_r;
%     model2ROIstruct_mfx_sep{i,1}.roiName=maskName;
%             userOptions.figureIndex = [4*(i-1)+3 (4*(i-1))+4];
% 
%     userOptions.figure1filename = [userOptions.analysisName '_' roiName '_P_barGraph'];
%     userOptions.figure2filename = [userOptions.analysisName '_' roiName '_P_Pvals'];
%     stats_p_r=compareRefRDM2candRDMs(pRDMs{i}, pSmodels, userOptions);
%     close all;
%     model2ROIstruct_mfx_sep{i,2}=stats_p_r;
%     model2ROIstruct_mfx_sep{i,2}.roiName=roiName;
% end
% 
% save([userOptions.analysisName '_ROI_ISO.mat'], '-v7.3')

end


