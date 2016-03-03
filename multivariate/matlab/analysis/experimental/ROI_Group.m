function [] = ROI_Group()

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

try
    toolboxRoot = 'D:\GitHub\LanguageMVPA\multivariate\matlab'; addpath(genpath(toolboxRoot));
catch
    toolboxRoot = '~/GitHub/LanguageMVPA/multivariate/matlab'; addpath(genpath(toolboxRoot));
end
% cd /space/raid5/data/monti/Analysis/LanguageMVPA/RSA
cd D:\fmri\LanguageMVPA
% userOptions = defineUserOptions_Group();
userOptions = defineUserOptions_LSA;
gotoDir(userOptions.rootPath);
% userOptions.analysisType=analysisType;
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
% %%%%%%%%%%%%%%%%%%%%%%
% % load data and masks
% fullBrainVols=fMRIDataPreparation(betaCorrespondence_LSA(), userOptions);
% binaryMasks_nS = fMRIMaskPreparation(userOptions);
% responsePatterns = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence_LSA(), userOptions);
% clear fullBrainVols binaryMasks_nS


userOptions.analysisName='tstat';

load(['ImageData/' userOptions.analysisName '_responsePatterns.mat'])

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
% %%%%%%%%%%%%%%%%%%%%%%
syntax=[repmat([4 3 2 1], 1, 8) repmat(4+[4 3 2 1], 1, 8)];
verb=[ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)];
verb = [verb 8+verb];
anim = [ones(1,16) 2*ones(1,16) 3*ones(1,16) 4*ones(1,16)];
% take mean of each condition to visualize...

semanticsResponsePatterns=responsePatterns;
syntaxResponsePatterns=responsePatterns;
animResponsePatterns=responsePatterns;
for m = 1:length(userOptions.maskNames)
    currMask=userOptions.maskNames{m};
    for s = 1:length(userOptions.subjectNames)
        
        currSub=userOptions.subjectNames{s};
        currPat=responsePatterns.(currMask).(currSub);
        semanticsMat=zeros(max(verb),size(currPat,1));
        syntaxMat=zeros(max(syntax),size(currPat,1));
        animMat=zeros(max(anim),size(currPat,1));
        for i = 1:max(verb)
            semanticsMat(i,:) = mean(currPat(:, verb==i)');
        end
        
        for i = 1:max(syntax)
            syntaxMat(i,:) = mean(currPat(:, syntax==i)');
        end
        
        for i = 1:max(anim)
            animMat(i,:) = mean(currPat(:, anim==i)');
        end
        semanticsResponsePatterns.(currMask).(currSub) = semanticsMat';
        syntaxResponsePatterns.(currMask).(currSub) = syntaxMat';
        animResponsePatterns.(currMask).(currSub) = animMat';
        %         semanticsResponsePatterns.(currMask).(currSub)=currPat(:,1:16);
        %         sSemanticsResponsePatterns.(currMask).(currSub)=currPat(:,1:8);
        %         pSemanticsResponsePatterns.(currMask).(currSub)=currPat(:,9:16);
        %         syntaxResponsePatterns.(currMask).(currSub)=currPat(:,17:24);
        %         sSyntaxResponsePatterns.(currMask).(currSub)=currPat(:,17:20);
        %         pSyntaxResponsePatterns.(currMask).(currSub)=currPat(:,21:24);
        %         animResponsePatterns.(currMask).(currSub)=currPat(:,25:28);
        %         sAnimResponsePatterns.(currMask).(currSub)=currPat(:,25:26);
        %         pInanimResponsePatterns.(currMask).(currSub)=currPat(:,27:28);
    end
end

%
% %% RSA %%
% %%%%%%%%%%%%%%%%%%%%%
% %% RDM calculation %%
% %%%%%%%%%%%%%%%%%%%%%
disp(['Starting RSA analysis']);
condLabels_bac = userOptions.conditionLabels;


userOptions.analysisName='avg_semantics';
ssemanticsRDMs = constructRDMs(semanticsResponsePatterns, betaCorrespondence_semantics(), userOptions);
semanticsRDMs = averageRDMs_subjectSession(ssemanticsRDMs, 'session', 'subject');
userOptions.conditionLabels={'s_Touch', 's_Light', 's_Hit', 's_Crush', ...
    's_Kiss', 's_Stretch', 's_Kick', 's_Console', 'l_Touch', ...
    'l_Light', 'l_Hit', 'l_Crush', 'l_Kiss', 'l_Stretch', 'l_Kick', 'l_Console' };
MDSConditions(semanticsRDMs, userOptions);

userOptions.analysisName='avg_syntax';
ssyntaxRDMs = constructRDMs(syntaxResponsePatterns, betaCorrespondence_syntax(), userOptions);
syntaxRDMs = averageRDMs_subjectSession(ssyntaxRDMs, 'session', 'subject');
userOptions.conditionLabels = {'S_PR', 'S_PC', 'S_AR', 'S_AC', 'P_PR', 'P_PC', 'P_AR', 'P_AC'};
MDSConditions(syntaxRDMs, userOptions);

userOptions.analysisName='avg_anim';
sanimRDMs = constructRDMs(animResponsePatterns, betaCorrespondence_anim(), userOptions);
animRDMs = averageRDMs_subjectSession(sanimRDMs, 'session', 'subject');
userOptions.conditionLabels = {'S_Inanim', 'S_anim', 'P_inanim', 'P_anim'};
MDSConditions(animRDMs, userOptions);


figureRDMs(semanticsRDMs, userOptions, struct('fileName', 'Semantics_RoIRDMs', 'figureNumber', 1));
figureRDMs(syntaxRDMs, userOptions, struct('fileName', 'Syntax_RoIRDMs', 'figureNumber', 2));
figureRDMs(animRDMs, userOptions, struct('fileName', 'Anim_RoIRDMs', 'figureNumber', 3));

% load(['RDMs/' userOptions.analysisName '_RDMs.mat']);

% sRDMs = averageRDMs_subjectSession(RDMs, 'session');
% RDMs = averageRDMs_subjectSession(RDMs, 'session', 'subject');
semModels = constructModelRDMs(modelRDMs_avg_semantics(),userOptions);
synModels = constructModelRDMs(modelRDMs_avg_syntax(),userOptions);

semLModels = semModels;
semPModels = semModels;

for i = 1:length(semModels)
    semLModels(i).RDM(9:16, 9:16) = nan;
%     semLModels(i).RDM(:, 9:16) = nan;
    semPModels(i).RDM(1:8, 1:8) = nan;
%     semPModels(i).RDM(1:8, :) = nan;
    semLModels(i).RDM(logical(eye(length(semLModels(i).RDM)))) = 0;
    semPModels(i).RDM(logical(eye(length(semLModels(i).RDM)))) = 0;
end
synLModels = synModels;
synPModels = synModels;

for i = 1:length(synModels)
    synLModels(i).RDM(5:8, 5:8) = nan;
%     synLModels(i).RDM(:, 5:8) = nan;
    synPModels(i).RDM(1:4, 1:4) = nan;
%     synPModels(i).RDM(1:4, :) = nan;
    synLModels(i).RDM(logical(eye(length(synLModels(i).RDM)))) = 0;
    synPModels(i).RDM(logical(eye(length(synLModels(i).RDM)))) = 0;
end

figureRDMs(semModels, userOptions, struct('fileName', 'Semantics_ModelRDMs', 'figureNumber', 4));
figureRDMs(synModels, userOptions, struct('fileName', 'Syntax_ModelRDMs', 'figureNumber', 5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dendrogramConditions(sSemanticsRDMs, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% relationship amongst multiple RDMs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pairwiseCorrelateRDMs({RDMs, Models}, userOptions);
% pairwiseCorrelateRDMs({sentRDMs, sentModels}, userOptions);
MDSRDMs({semanticsRDMs, semModels}, userOptions);
MDSRDMs({semanticsRDMs, semLModels}, userOptions);
MDSRDMs({semanticsRDMs, semPModels}, userOptions);
MDSRDMs({syntaxRDMs, synModels}, userOptions);
MDSRDMs({syntaxRDMs, synLModels}, userOptions);
MDSRDMs({syntaxRDMs, synPModels}, userOptions);
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
for i = 1: size(ssemanticsRDMs,1)
    asemanticsRDMs{i}=ssemanticsRDMs(i,:);
    asyntaxRDMs{i}=ssyntaxRDMs(i,:);
end
for i=1:numel(semModels)
    lSemmodels{i}=semLModels(i);
    pSemmodels{i}=semPModels(i);
end
for i=1:numel(synModels)
    lSynmodels{i}=synLModels(i);
    pSynmodels{i}=synPModels(i);
end
% for i=1:numel(Models)
%     models{i}=Models(i);
% end
clear i
close all

for i = 1:numel(asemanticsRDMs)
    userOptions.figureIndex = [4*(i-1)+1 (4*(i-1))+2];
    
    roiName=userOptions.maskNames{i};
    disp(['Processing ' roiName])
    %     %set other mode model RDM to empirical
    userOptions.figure1filename = [userOptions.analysisName '_semantics_' roiName '_barGraph'];
    userOptions.figure2filename = [userOptions.analysisName '_semantics_' roiName '_Pvals'];
    %     %         switched so we match ROIs to model
    %     %             stats_p_r.(modelName)=compareRefRDM2candRDMs(models, RDMs(roiIndex), userOptions);
    %     userOptions.resultsPath = fullfile(userOptions.rootPath, 'Statistics/');
    stats_p_r=compareRefRDM2candRDMs(asemanticsRDMs{i}, lSemmodels, userOptions);
    semantics_model2ROIstruct_mfx{i,1}=stats_p_r;
    semantics_model2ROIstruct_mfx{i,1}.roiName=roiName;
    stats_p_r=compareRefRDM2candRDMs(asemanticsRDMs{i}, pSemmodels, userOptions);
    semantics_model2ROIstruct_mfx{i,2}=stats_p_r;
    semantics_model2ROIstruct_mfx{i,2}.roiName=roiName;
    
    
    userOptions.figureIndex = [4*(i-1)+3 (4*(i-1))+4];

    stats_p_r=compareRefRDM2candRDMs(asyntaxRDMs{i}, lSynmodels, userOptions);
    syntax_model2ROIstruct_mfx{i,1}=stats_p_r;
    syntax_model2ROIstruct_mfx{i,1}.roiName=roiName;
    stats_p_r=compareRefRDM2candRDMs(asyntaxRDMs{i}, pSynmodels, userOptions);
    syntax_model2ROIstruct_mfx{i,2}=stats_p_r;
    syntax_model2ROIstruct_mfx{i,2}.roiName=roiName;
end


save([userOptions.rootPath '/Statistics/avg_syntax_semantics_ROI.mat'], '-v7.3')
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


