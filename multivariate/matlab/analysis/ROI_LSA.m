function [] = ROI_LSA(analysisType)

% LMVPA_fMRI... depending on good beta extraction. This script should work
% I think...
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
    toolboxRoot = 'D:\GitHub\LanguageMVPA\multivariate\matlab'; addpath(genpath(toolboxRoot));
catch
    toolboxRoot = '~/GitHub/LanguageMVPA/multivariate/matlab'; addpath(genpath(toolboxRoot));
end
% cd /space/raid5/data/monti/Analysis/LanguageMVPA/RSA
cd D:\fmri\LanguageMVPA
% userOptions = defineUserOptions_Group();
userOptions = defineUserOptions_LSA;
gotoDir(userOptions.rootPath);
userOptions.analysisType=analysisType;
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

if strcmpi(userOptions.analysisType, 'SVM')
    
    % %% SVM %%
    % need a file that establishes runs and labels:
    syntax=repmat([4 3 2 1], 1, 16)';
    actpass=repmat([2 2 1 1], 1, 16)';
    relcan=repmat([2 1 2 1], 1, 16)';
    verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,2)';
    verbRuns=syntax;
    syntaxRuns=verb;
    % clear responsePatterns
    nPerm=100;
    theseRuns=syntaxRuns(1:32);
    thisLabel=permuteLabels(syntax(1:32), theseRuns, nPerm);
    OPTS=['-q -s 0 -t 2'];
    res=repmat(struct('acc', zeros(length(unique(theseRuns)), nPerm), ...
        'mse', zeros(length(unique(theseRuns)), nPerm), ...
        'scc', zeros(length(unique(theseRuns)), nPerm), ...
        'p95threshold', inf, 'classAcc', 0), ...
        length(userOptions.maskNames), length(userOptions.subjectNames));
    for iMask=1:length(userOptions.maskNames)
        disp(['Running ' userOptions.maskNames{iMask}]);
        for iSub = 1:length(userOptions.subjectNames)
            res(iMask, iSub) = CVpermutationSVM( ...
                responsePatterns.(userOptions.maskNames{iMask}).(userOptions.subjectNames{iSub})(:,1:32)', ...
                thisLabel, theseRuns, 'opts', OPTS);
        end
    end
    save(['Statistics/' userOptions.analysisName '_SVM.mat'], res, userOptions, OPTS);
elseif strcmpi(userOptions.analysisType, 'RSA')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %% RSA %%
    % %%%%%%%%%%%%%%%%%%%%%
    % %% RDM calculation %%
    % %%%%%%%%%%%%%%%%%%%%%
    % disp(['Starting RSA analysis']);
    try
        load(['RDMs/' userOptions.analysisName '_RDMs.mat']);
    catch
        RDMs = constructRDMs(responsePatterns, betaCorrespondence_LSA(), userOptions);
    end
    sRDMs = averageRDMs_subjectSession(RDMs, 'session');
    RDMs = averageRDMs_subjectSession(RDMs, 'session', 'subject');
    
    Models=constructModelRDMs(modelRDMs_Compute(),userOptions);
    LModels = Models;
    PModels = Models;
    % Language Models Only
    for i = 1:length(Models)
        LModels(i).RDM(33:64, :) = nan;
        LModels(i).RDM(:, 33:64) = nan;
        PModels(i).RDM(:, 1:32) = nan;
        PModels(i).RDM(1:32, :) = nan;
        LModels(i).RDM(logical(eye(length(LModels(i).RDM)))) = 0;
        PModels(i).RDM(logical(eye(length(PModels(i).RDM)))) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% First-order visualisation %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figureRDMs(RDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1));
    figureRDMs(LModels, userOptions, struct('fileName', 'ComputeRDMs', 'figureNumber', 2));
    figureRDMs(PModels, userOptions, struct('fileName', 'ComputeRDMs', 'figureNumber', 3));
    
    %     MDSConditions(RDMs, userOptions);
    %     dendrogramConditions(sentRDMs, userOptions);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% relationship amongst multiple RDMs %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % pairwiseCorrelateRDMs({RDMs, Models}, userOptions);
    % pairwiseCorrelateRDMs({sentRDMs, sentModels}, userOptions);
    MDSRDMs({RDMs, LModels(1:end-2)}, userOptions);
    MDSRDMs({RDMs, PModels(2:end-2)}, userOptions);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% statistical inference %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % userOptions.RDMcorrelationType='Kendall_taua';
    userOptions.RDMcorrelationType='Spearman';
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
    % for i=1:numel(LModels)
    %     lFmodels{i}=LModels(i);
    %     pFmodels{i}=PModels(i);
    % end
    for i=1:numel(Models)
        lmodels{i}=LModels(i);
        pmodels{i}=PModels(i);
    end
    clear i
    close all
    %%% separating L and P
    for i = 1:numel(aRDMs)
        userOptions.figureIndex = [4*(i-1)+1 (4*(i-1))+2];
        roiName=userOptions.maskNames{i};
        disp(['Processing ' roiName])
        %set other mode model RDM to empirical
        lmodels{end-1}.RDM(1:32,1:32) = RDMs(i).RDM(1:32,1:32);
        lmodels{end}.RDM(1:32, 1:32) = RDMs(i).RDM(33:64,33:64);
        pmodels{end-1}.RDM(33:64,33:64) = RDMs(i).RDM(1:32,1:32);
        pmodels{end}.RDM(33:64, 33:64) = RDMs(i).RDM(33:64,33:64);
        userOptions.figure1filename = [userOptions.analysisName '_' roiName '_L_barGraph'];
        userOptions.figure2filename = [userOptions.analysisName '_' roiName '_L_Pvals'];
        %         switched so we match ROIs to model
        %             stats_p_r.(modelName)=compareRefRDM2candRDMs(models, RDMs(roiIndex), userOptions);
        userOptions.resultsPath = fullfile(userOptions.rootPath, 'Statistics/');
        stats_p_r=compareRefRDM2candRDMs(aRDMs{i}, lmodels, userOptions);
        model2ROIstruct_mfx{i,1}=stats_p_r;
        model2ROIstruct_mfx{i,1}.roiName=roiName;
        userOptions.figureIndex = [4*(i-1)+3 (4*(i-1))+4];
        userOptions.figure1filename = [userOptions.analysisName '_' roiName '_P_barGraph'];
        userOptions.figure2filename = [userOptions.analysisName '_' roiName '_P_Pvals'];
        stats_p_r=compareRefRDM2candRDMs(aRDMs{i}, pmodels, userOptions);
        model2ROIstruct_mfx{i,2}=stats_p_r;
        model2ROIstruct_mfx{i,2}.roiName=roiName;
    end
elseif strcmpi(userOptions.analysisType, 'power')
    % derpy power analysis
    syntax=repmat([4 3 2 1], 1, 16)';
    verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,2)';
    verbRuns=syntax;
    syntaxRuns=verb;
    % clear responsePatterns
    nPerm=100;
    theseRuns=verbRuns(1:32);
    thisLabel=permuteLabels(verb(1:32), theseRuns, nPerm);
    randMask=userOptions.maskNames{randi(length(userOptions.maskNames))};
    randSub=userOptions.subjectNames{randi(length(userOptions.subjectNames))};
    trueData=responsePatterns.(randMask).(randSub)(:,1:32);
    i=0;
    powAnalysis=figure;
    for sigma=1:5:100
        i=i+1;
        randData = bsxfun(@plus, thisLabel(:,1), normrnd(0, sigma, size(trueData,2), size(trueData,1)));
        t = CVpermutationSVM(randData,thisLabel, theseRuns);
        subplot(4,5,i)
        histogram(mean(t.acc))
        line([t.classAcc t.classAcc], [0 30], 'Marker', 'o')
        line([t.p95threshold t.p95threshold], [0 30], 'Marker', '*')
        legend('Empirical dist', 'accuracy', 'threshold');
        title(['SNR: 1/ ' num2str(sigma)])
        [r(i),p(i)] = corr(randData(:), trueData(:));
    end
    figure;
    plot(1:5:100, r, 'o');
    title('Correlation of true and random data')
    xlabel('SNR')
    ylabel('Corr');
    return
else
    disp('Invalid analysis string')
    return
end
save([userOptions.rootPath '/Statistics/' userOptions.analysisName '_' userOptions.analysisType '_ROI.mat'], '-v7.3')


end
