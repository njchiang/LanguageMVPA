function [] = ROI_Music(analysisType)

% Music_fMRI
% this script performs region of interest analysis on fMRI data.
% Written by Jeff Chiang
% Based on Cai Wingfield 5-2010, 6-2010, 7-2010, 8-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

toolboxRoot = '/space/raid5/data/monti/Analysis/LanguageMVPA/Music/mvpa_data/code'; addpath(genpath(toolboxRoot));
userOptions = defineUserOptions_music();
gotoDir(userOptions.rootPath);
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%
% load data and masks
% fullBrainVols = fMRIDataPreparation('SPM', userOptions);
fullBrainVols=fMRIDataPreparation(betaCorrespondence_music(), userOptions);
binaryMasks_nS = fMRIMaskPreparation(userOptions);
responsePatterns = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence_music(), userOptions);
clear fullBrainVols binaryMasks_nS

%% SVM %%
if strcmp(analysisType, 'SVM')
    %initialize
    outDir=fullfile(userOptions.rootPath, 'Statistics', 'SVM');
    opts=['-s 1 -v 2 -t 0 -q']; %split data in half, linear SVM, etc.
    topts=['-s 1 -t 0 -q'];
    
    models=makeLabels_music();
    LMStimType=models(1).label;
    LMSyntax=models(2).label;
    LSyntax=LMSyntax(LMStimType==1);
    MSyntax=LMSyntax(LMStimType==2);
    ActPass=models(6).label;
    LMStimType_acc=zeros(length(userOptions.subjectNames), length(userOptions.maskNames));
    LMSyntax_acc=zeros(length(userOptions.subjectNames), length(userOptions.maskNames));
    LSyntax_acc=zeros(length(userOptions.subjectNames), length(userOptions.maskNames));
    MSyntax_acc=zeros(length(userOptions.subjectNames), length(userOptions.maskNames));
    LActPass_acc=zeros(length(userOptions.subjectNames), length(userOptions.maskNames));
    MActPass_acc=zeros(length(userOptions.subjectNames), length(userOptions.maskNames));
    LvMSyntax_acc=zeros(length(userOptions.subjectNames), length(userOptions.maskNames));
    MvLSyntax_acc=zeros(length(userOptions.subjectNames), length(userOptions.maskNames));
    masks=fieldnames(responsePatterns);
    for m = 1:length(masks)
        currMask=masks{m};
        currMaskData=responsePatterns.(currMask);
        currSubjs=fieldnames(currMaskData);
        for s = 1:length(currSubjs)
            currSub=currSubjs{s};
            disp(['Processing ' currSub ' ' currMask ]);
            currPattern=currMaskData.(currSub);
            currPattern=zscore(currPattern,0,2);
            L_patterns=currPattern(:,LMStimType==1);
            M_patterns=currPattern(:,LMStimType==2);
            LMStimType_acc(s, m)=libsvmtrain(LMStimType',currPattern',opts);
            LMSyntax_acc(s,m) = libsvmtrain(LMSyntax', currPattern', opts);
            LSyntax_acc(s,m) = libsvmtrain(LMSyntax(LMStimType==1)', L_patterns', opts);
            MSyntax_acc(s,m)=libsvmtrain(LMSyntax(LMStimType==2)', M_patterns', opts);
            LActPass_acc(s,m)=libsvmtrain(ActPass(ActPass~=0)', L_patterns(ActPass~=0)',opts);
            MActPass_acc(s,m)=libsvmtrain(ActPass(ActPass~=0)', M_patterns(ActPass~=0)',opts);
            LvMStruct=libsvmtrain(LMSyntax(LMStimType==1)', L_patterns',topts);
            [~, acc, ~] = libsvmpredict(LMSyntax(LMStimType==2)',M_patterns', LvMStruct);
            LvMSyntax_acc(s,m)=acc(1);
            MvLStruct=libsvmtrain(LMSyntax(LMStimType==2)', M_patterns',topts);
            [~, acc, ~] = libsvmpredict(LMSyntax(LMStimType==1)',L_patterns', MvLStruct);
            MvLSyntax_acc(s,m)=acc(1);
        end
    end
    SVM_ROI_results.StimType=LMStimType_acc;
    SVM_ROI_results.Cross_Stimulus_Syntax=LMSyntax_acc;
    SVM_ROI_results.Language_Syntax=LSyntax_acc;
    SVM_ROI_results.Music_Syntax=MSyntax_acc;
    SVM_ROI_results.L_to_M_Syntax=LvMSyntax_acc;
    SVM_ROI_results.M_to_L_Syntax=MvLSyntax_acc;
    SVM_ROI_results.LActPass=LActPass_acc;
    SVM_ROI_results.MActPass=MActPass_acc;
    returnHere=pwd;
    gotoDir(outDir)
    csvwrite('StimType.csv', LMStimType_acc)
    csvwrite('Cross_Stimulus_Syntax.csv', LMSyntax_acc);
    csvwrite('Language_Syntax.csv', LSyntax_acc);
    csvwrite('Music_Syntax.csv', MSyntax_acc);
    csvwrite('L_to_M_Syntax.csv', LvMSyntax_acc);
    csvwrite('M_to_L_Syntax.csv', MvLSyntax_acc);
    csvwrite('L_ActvPass.csv', LActPass_acc);    
    csvwrite('M_ActvPass.csv', MActPass_acc);    
    save(['SVM_ROI_Results.mat'], 'SVM_ROI_results');
    gotoDir(returnHere);
    
    %% RSA %%
elseif strcmp(analysisType, 'RSA')
    %%%%%%%%%%%%%%%%%%%%%
    %% RDM calculation %%
    %%%%%%%%%%%%%%%%%%%%%
    
    RDMs = constructRDMs(responsePatterns, betaCorrespondence_music(), userOptions);
    sRDMs = averageRDMs_subjectSession(RDMs, 'session');
    RDMs = averageRDMs_subjectSession(RDMs, 'session', 'subject');
    
    Models = constructModelRDMs(modelRDMs_music(), userOptions);
    
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
    %     MDSRDMs({RDMs, Models}, userOptions);
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
    userOptions.candRDMdifferencesMultipleTesting = 'none';
    
    %probably need to average...
    % roiIndex = 1;% index of the ROI for which the group average RDM will serve
    % as the reference RDM.
    %need to reform RDMs as cell array
    for i = 1: numel(RDMs)
        aRDMs{i}=RDMs(i);
    end
    for i=1:numel(Models)
        models{i}=Models(i);
    end
    clear i
    for i=1:numel(Models)
        modelName=Models(i).name;
        disp(['Processing ' modelName])
        userOptions.resultsPath = fullfile(userOptions.rootPath, 'Statistics');
        userOptions.figure1filename = [modelName '_FDR__barGraph'];
        userOptions.figure2filename = [modelName '_FDR_Pvals'];
        % switched so we match ROIs to model
        %     stats_p_r.(modelName)=compareRefRDM2candRDMs(models, RDMs(roiIndex), userOptions);
        stats_p_r.(modelName)=compareRefRDM2candRDMs(models{i}, aRDMs, userOptions);
        close all
    end
    
    save([userOptions.rootPath '/Statistics/RSA_ROI_analysis.mat'], '-v7.3')
else
    disp('Not a valid classifier, exiting...')
end
end
