
% toolboxRoot = '/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/code'; 
% addpath(genpath(toolboxRoot));
userOptions = defineUserOptions_behavioral();
% gotoDir(userOptions.rootPath);
analysisType='RSA';
userOptions.analysisType=analysisType;
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%

% load subject.mat
% mat = 1-results;
mat = rand(32,32);
mat = mat+mat';
cond = 'test';
subID = '001';
%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

% RDMs = constructRDMs(responsePatterns, betaCorrespondence_RSA(), userOptions);
RDMs(1,1) = rdmWrapper(mat, cond, subID);
sRDMs = averageRDMs_subjectSession(RDMs, 'session');
RDMs = averageRDMs_subjectSession(RDMs, 'session', 'subject');

Models=constructModelRDMs(modelRDMs_behavioral(),userOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figureRDMs(RDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1));
figureRDMs(Models, userOptions, struct('fileName', 'ModelRDMs', 'figureNumber', 2));

%     MDSConditions(RDMs, userOptions);
%     dendrogramConditions(RDMs, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% relationship amongst multiple RDMs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pairwiseCorrelateRDMs({RDMs, Models}, userOptions);
MDSRDMs({RDMs, Models}, userOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistical inference %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

userOptions.RDMcorrelationType='Kendall_taua';
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
userOptions.plotpValues='*';


%probably need to average...
% roiIndex = 1;% index of the ROI for which the group average RDM will serve
% as the reference RDM.
%need to reform RDMs as cell array
for i = 1: size(sRDMs,1)
    aRDMs{i}=sRDMs(i,:);
end
for i=1:numel(Models)
    models{i}=Models(i);
end
clear i
close all

for i = 1:numel(aRDMs)
    roiName=RDMs(i).name;
    disp(['Processing ' roiName])
    userOptions.resultsPath = fullfile(userOptions.rootPath, 'Statistics');
    userOptions.figure1filename = [roiName '_FDR_mfx_barGraph'];
    userOptions.figure2filename = [roiName '_FDR_mfx_Pvals'];
    userOptions.RDMcorrelationType='Kendall_taua';
    %         switched so we match ROIs to model
    %             stats_p_r.(modelName)=compareRefRDM2candRDMs(models, RDMs(roiIndex), userOptions);
    stats_p_r{i}=compareRefRDM2candRDMs(aRDMs{i}, models, userOptions);
    close all;
end
