function [] = ROI_259()

% LMVPA_fMRI
%preprocessing was all done, so just go straight to loading RDMs!
% all comparisons will be done on sentence stimuli ONLY
% this script performs region of interest analysis on fMRI data.
% Written by Jeff Chiang
% Based on Cai Wingfield 5-2010, 6-2010, 7-2010, 8-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council
% comparing within vs between:

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

% toolboxRoot = '/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/code'; addpath(genpath(toolboxRoot));
% cd /space/raid5/data/monti/Analysis/LanguageMVPA/RSA
% userOptions = defineUserOptions_ROI();
userOptions = defineUserOptions_259();
gotoDir(userOptions.rootPath);
analysisType='RSA';

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
% %%%%%%%%%%%%%%%%%%%%%%
% % load data and masks
% % fullBrainVols = fMRIDataPreparation('SPM', userOptions);
% fullBrainVols=fMRIDataPreparation(betaCorrespondence_LMVPA(), userOptions);
% binaryMasks_nS = fMRIMaskPreparation(userOptions);
% responsePatterns = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence_LMVPA(), userOptions);
% 
% clear fullBrainVols binaryMasks_nS
%
% %% RSA %%
% %%%%%%%%%%%%%%%%%%%%%
% %% RDM calculation %%
% %%%%%%%%%%%%%%%%%%%%%
disp(['Starting RSA analysis']);
% RDMs = constructRDMs(responsePatterns, betaCorrespondence_LMVPA(), userOptions);

load('../Topics/259_RDMs.mat')


sRDMs = averageRDMs_subjectSession(RDMs, 'session');
RDMs = averageRDMs_subjectSession(RDMs, 'session', 'subject');

%sentences only:
RDMs.RDM=RDMs.RDM(1:32, 1:32);
% Models=constructModelRDMs(modelRDMs_ROI(),userOptions);
% Models=constructModelRDMs(modelRDMs_ROI_edit(), userOptions);
% Models=constructModelRDMs(modelRDMs_ROI_layered(), userOptions);
Models=constructModelRDMs(modelRDMs_259(), userOptions);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figureRDMs(RDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1));
figureRDMs(Models, userOptions, struct('fileName', 'ModelRDMs', 'figureNumber', 2));

%     MDSConditions(sentRDMs, userOptions);
%     dendrogramConditions(sentRDMs, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% relationship amongst multiple RDMs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pairwiseCorrelateRDMs({RDMs, Models}, userOptions);

%     MDSRDMs({RDMs, Models}, userOptions);
%     MDSRDMs({sentRDMs, sentModels}, userOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistical inference %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% userOptions.RDMcorrelationType='Kendall_taua';
userOptions.RDMcorrelationType='Spearman';

userOptions.RDMrelatednessTest = 'subjectRFXsignedRank';
% userOptions.RDMrelatednessTest = 'conditionRFXbootstrap';

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
userOptions.barsOrderedByRDMCorr='false';


for i=1:numel(Models)
    models{i}=Models(i);
end
% stats_p_r=compareRefRDM2candRDMs(subjectRDMs, modelRDMs_cell, userOptions);


    userOptions.figure1filename = ['TopicsProjectBarGraph'];
    userOptions.figure2filename = [ 'TopicsProjectPVals'];
    stats_p_r=compareRefRDM2candRDMs(RDMs, models, userOptions);
    close all
end


save([userOptions.resultsPath '/' userOptions.analysisName '_Model2ROI_bootstrap.mat'], '-v7.3')

end
