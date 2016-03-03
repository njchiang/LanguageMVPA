function [] = ROI_LSA_SVM()

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

% toolboxRoot = '~/Box Sync/UCLA/Research/LanguageMVPA/code'; addpath(genpath(toolboxRoot));
try
    toolboxRoot = 'Z:\GitHub\LanguageMVPA\multivariate\matlab'; addpath(genpath(toolboxRoot));
catch
    toolboxRoot = '~/GitHub/LanguageMVPA/multivariate/matlab'; addpath(genpath(toolboxRoot));
end
% cd /space/raid5/data/monti/Analysis/LanguageMVPA/RSA
cd Z:\fmri\LanguageMVPA
% userOptions = defineUserOptions_Group();
userOptions = defineUserOptions_LSA;
gotoDir(userOptions.rootPath);
analysisType='SVM';
userOptions.analysisType=analysisType;
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
% %%%%%%%%%%%%%%%%%%%%%%
% % load data and masks
% fullBrainVols=fMRIDataPreparation(betaCorrespondence_LSA(), userOptions);
% binaryMasks_nS = fMRIMaskPreparation(userOptions);
% responsePatterns = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence_LSA(), userOptions);
% clear fullBrainVols binaryMasks_nS

load('ImageData/LSA_responsePatterns.mat')
% for generic
% masks = fieldnames(responsePatterns);
% subjs = fieldnames(responsePatterns.(masks{1}));
% need a file that establishes runs and labels:
syntax=repmat([4 3 2 1], 1, 16)';
actpass=repmat([2 2 1 1], 1, 16)';
relcan=repmat([2 1 2 1], 1, 16)';
verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,2)';

verbRuns=syntax;
syntaxRuns=verb;
% clear responsePatterns
thisLabel=syntax(1:32);
theseRuns=syntaxRuns(1:32);
OPTS=['-q -s 0 -t 2'];
CVOPTS=[OPTS ' -v 4'];
nPerm=1;
fullZr=zeros(length(userOptions.subjectNames),1);
fullZp=fullZr;
chunkZr=fullZr;
chunkZp=fullZr;
res=repmat(struct('acc', zeros(length(unique(theseRuns)), nPerm), ...
    'mse', zeros(length(unique(theseRuns)), nPerm), ...
    'scc', zeros(length(unique(theseRuns)), nPerm), ...
    'p95threshold', inf, 'classAcc', 0), ...
    length(userOptions.maskNames), length(userOptions.subjectNames));
for iMask=1:length(userOptions.maskNames)
    for iSub = 1:length(userOptions.subjectNames)
        res(iSub) = CVpermutationSVM( ...
            responsePatterns.(mask).(userOptions.subjectNames{iSub})(:,1:32)', ...
            thisLabel, theseRuns, 'opts', OPTS);
    end
end
save([userOptions.rootPath '/Statistics/' userOptions.analysisName '_ROI.mat'], '-v7.3')


end
%         % check correlation of full zscore
%         [fullZr(iSub), fullZp(iSub)] =corr(pdist(corr(zscore( ...
%             responsePatterns.(mask).(userOptions.subjectNames{iSub})(:, 1:32),0,2)))',...
%             pdist(corr( ...
%             responsePatterns.(mask).(userOptions.subjectNames{iSub})(:, 1:32)))');
%         % check correlation of our zscore
%         [chunkZr(iSub), chunkZp(iSub)] =corr(pdist(corr( ...
%             responsePatterns.(mask).(userOptions.subjectNames{iSub})(:,1:32)))',...
%             pdist(corr( ...
%             responsePatterns.(mask).(userOptions.subjectNames{iSub})(:,1:32)))');
