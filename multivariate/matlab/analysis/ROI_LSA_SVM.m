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

toolboxRoot = 'Z:/Box Sync/UCLA/Research/LanguageMVPA/code'; addpath(genpath(toolboxRoot));
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

a=load(['ImageData/ROI_LSA_responsePatterns.mat']);
responsePatterns=a.responsePatterns;
clear a
% don't need to separate before correlation, they will be the same rank
% order after separating
models=makeLabels_LSA;
% models of interest: 2 3 4 7 8 9
% runs=ones(32,1);
nPerm=1000;
svmKernel=0;
verbRuns=models(10).label(1:32);
syntaxRuns=models(2).label(1:32);
% syntax-- within
syntaxLabel=models(10).label(1:32);
syntaxPermutedLabels=permuteLabels(syntaxLabel, syntaxRuns, nPerm);
actpassLabel=models(11).label(1:32);
actpassPermutedLabels=permuteLabels(actpassLabel, syntaxRuns, nPerm);
relcanLabel=models(12).label(1:32);
relcanPermutedLabels=permuteLabels(relcanLabel, syntaxRuns, nPerm);
arpcLabel=models(13).label(1:32);
arpcPermutedLabels=permuteLabels(arpcLabel(arpcLabel>0), syntaxRuns(arpcLabel>0), nPerm);

verbLabel=models(2).label(1:32);
verbPermutedLabels=permuteLabels(verbLabel, verbRuns, nPerm);
animFamLabel=models(4).label(1:32);
animFamPermutedLabels=permuteLabels(animFamLabel(animFamLabel>0), verbRuns(animFamLabel>0), nPerm);
inAnimFamLabel=models(5).label(1:32);
inAnimFamPermutedLabels=permuteLabels(inAnimFamLabel(inAnimFamLabel>0), verbRuns(inAnimFamLabel>0), nPerm);
animLabel=models(3).label(1:32);
animPermutedLabels=permuteLabels(animLabel, verbRuns, nPerm);
% 
% complexRuns=syntaxRuns(1:16);
% complexLabel=models(15).label(1:32);
% balComplexLabel=[2*ones(8,1); ones(8,1)];
% complexPermutedLabels=permuteLabels(balComplexLabel, syntaxRuns(1:16) ,nPerm);
h=waitbar(0, ['Processing ...']);

for m = 1:length(userOptions.maskNames)-1
    currMask=userOptions.maskNames{m};
    tl2l=[];
    tl2p=[];
    tp2p=[];
    tp2l=[];
    parfor s = 1:length(userOptions.subjectNames)
        currSub=userOptions.subjectNames{s};
%         currPat=zscore(responsePatterns.(currMask).(currSub))';
%         currPat=zscore(responsePatterns.(currMask).(currSub)');
        currPat=(responsePatterns.(currMask).(currSub)');
%         zCurrPat=zscore((responsePatterns.(currMask).(currSub)'));
        [tl2l(s).syntax, tl2p(s).syntax, tp2p(s).syntax, tp2l(s).syntax, lsyntax(s,m), psyntax(s,m)] ...
            = run_SVM_LSA_LMVPA(syntaxLabel, syntaxPermutedLabels, currPat, syntaxRuns, svmKernel);
        [tl2l(s).actpass, tl2p(s).actpass, tp2p(s).actpass, tp2l(s).actpass, lactpass(s,m), pactpass(s,m)] ...
            = run_SVM_LSA_LMVPA(actpassLabel, actpassPermutedLabels, currPat, syntaxRuns, svmKernel);
        [tl2l(s).relcan, tl2p(s).relcan, tp2p(s).relcan, tp2l(s).relcan, lrelcan(s,m), prelcan(s,m)] ...
            = run_SVM_LSA_LMVPA(relcanLabel, relcanPermutedLabels, currPat, syntaxRuns, svmKernel);
        [tl2l(s).verb, tl2p(s).verb, tp2p(s).verb, tp2l(s).verb, lverb(s,m), pverb(s,m)] ...
            = run_SVM_LSA_LMVPA(verbLabel, verbPermutedLabels, currPat, verbRuns, svmKernel);
%         [tl2l(s).fam, tl2p(s).fam, tp2p(s).fam, tp2l(s).fam, lfam(s,m), pfam(s,m)] ...
%             = run_SVM_LSA_LMVPA(famLabel, famPermutedLabels, currPat, verbRuns, svmKernel);
        [tl2l(s).anim, tl2p(s).anim, tp2p(s).anim, tp2l(s).anim, lanim(s,m), panim(s,m)] ...
            = run_SVM_LSA_LMVPA(animLabel, animPermutedLabels, currPat, verbRuns, svmKernel);
        [tl2l(s).animFam, tl2p(s).animFam, tp2p(s).animFam, tp2l(s).animFam, lanimFam(s,m), panimFam(s,m)] ...
            = run_SVM_LSA_LMVPA(animFamLabel, animFamPermutedLabels, currPat, verbRuns(animFamLabel>0), svmKernel);
        [tl2l(s).inAnimFam, tl2p(s).inAnimFam, tp2p(s).inAnimFam, tp2l(s).inAnimFam, linAnimFam(s,m), pinAnimFam(s,m)] ...
            = run_SVM_LSA_LMVPA(inAnimFamLabel, inAnimFamPermutedLabels, currPat, verbRuns(inAnimFamLabel>0), svmKernel);
        [tl2l(s).arpc, tl2p(s).arpc, tp2p(s).arpc, tp2l(s).arpc, larpc(s,m), parpc(s,m)] ...
            = run_SVM_LSA_LMVPA(arpcLabel, arpcPermutedLabels, currPat, syntaxRuns(arpcLabel>0), svmKernel);
    end
    l2l(:,m)=tl2l;
    l2p(:,m)=tl2p;
    p2p(:,m)=tp2p;
    p2l(:,m)=tp2l;
    waitbar(m/length(userOptions.maskNames));
end
close(h)
delete(gcp)
save('Statistics/Linear_zTime_ROI_SVM.mat', 'l2l', 'l2p', 'p2p', 'p2l', ...
    'lanim', 'panim', 'lanimFam', 'panimFam', ...
    'linAnimFam', 'pinAnimFam', 'lverb', 'pverb', ...
    'larpc', 'parpc', 'lsyntax', 'psyntax', 'lrelcan', 'prelcan', ...
    'userOptions')
% save('ROI_SVM_noComplex.mat', 'l2l', 'l2p', 'p2p', 'p2l', ...
%     'lanim', 'panim', 'lanimFam', 'panimFam', ...
%     'linAnimFam', 'pinAnimFam', 'lfam', 'pfam', 'lverb', 'pverb', ...
%     'larpc', 'parpc', 'lsyntax', 'psyntax', 'lrelcan', 'prelcan', ...
%     'userOptions')
% h=waitbar(0, ['Processing Complex...']);
% for m = 1:length(userOptions.maskNames)-1
%     currMask=userOptions.maskNames{m};
%     for s=1:length(userOptions.subjectNames)
%         currSub=userOptions.subjectNames{s};
%         currPat=zscore(responsePatterns.(currMask).(currSub))';
%         %complex has funky numbers, so do that individually... OR DOES THAT
%         %MATTER BECAUSE I DO IT NONPARMAETRICALLY?
%         sentPats=currPat(1:32, :);
%         picPats=currPat(33:64,:);
%         nR=10;
%         tmpl2lacc=[];
%         tmpl2ldist=[];
%         tmpl2lthresh=[];
%         tmpl2pacc=[];
%         tmpl2pdist=[];
%         tmpl2pthresh=[];
%         tmpp2pacc=[];
%         tmpp2pdist=[];
%         tmpp2pthresh=[];
%         tmpp2lacc=[];
%         tmpp2ldist=[];
%         tmpp2lthresh=[];
%         for i=1:nR
%             balanceSet=Shuffle(find(complexLabel==2))';
%             acSet=find(complexLabel==1)';
%             lSet=sentPats([balanceSet(1:8) acSet],:);
%             pSet=picPats([balanceSet(1:8) acSet],:);
%             [tmpl2lacc(i), tmpD, tmpl2lthresh(i)]=CVpermutationSVM(lSet, complexPermutedLabels, complexRuns);
%             tmpl2ldist(i,:)=mean(tmpD);
%             lComplex(s,m)=optimizeSVM(complexPermutedLabels(:,1), lSet, svmKernel);
%             if svmKernel==2
%                 opts=['-s 0 -t ' num2str(svmKernel) ' -c ' num2str(lComplex(s,m).best_C) ' -g ' num2str(lComplex(s,m).best_gamma)];
%             else
%                 opts=['-s 1 -t ' num2str(svmKernel) ' -n ' num2str(lComplex(s,m).best_C) ];
%             end
%             lStruct=libsvmtrain(complexPermutedLabels(:,1), lSet, opts);
%             [tmpl2pacc(i), tmpD, tmpl2pthresh(i)] = permutationSVM(lStruct, pSet, complexPermutedLabels, complexRuns);
%             tmpl2pdist(i,:) = mean(tmpD);
%             [tmpp2pacc(i), tmpD, tmpp2pthresh(i)]=CVpermutationSVM(pSet, complexPermutedLabels, complexRuns);
%             tmpp2pdist(i,:)=mean(tmpD);
%             pComplex(s,m)=optimizeSVM(complexPermutedLabels(:,1), pSet, svmKernel);
%             if svmKernel==2
%                 opts=['-s 0 -t ' num2str(svmKernel) ' -c ' num2str(pComplex(s,m).best_C) ' -g ' num2str(pComplex(s,m).best_gamma)];
%             else
%                 opts=['-s 1 -t ' num2str(svmKernel) ' -n ' num2str(pComplex(s,m).best_C) ];
%             end
%             pStruct=libsvmtrain(complexPermutedLabels(:,1), pSet, opts);
%             [tmpp2lacc(i), tmpD, tmpp2lthresh(i)] = permutationSVM(pStruct, lSet, complexPermutedLabels, complexRuns);
%             tmpp2ldist(i,:)=mean(tmpD);
%             tmpD=[];
%         end
%         
%         
%         l2l(s,m).Complex.acc=mean(tmpl2lacc);
%         l2p(s,m).Complex.acc=mean(tmpl2pacc);
%         p2p(s,m).Complex.acc=mean(tmpp2lacc);
%         p2l(s,m).Complex.acc=mean(tmpp2pacc);
%         l2l(s,m).Complex.dist=mean(tmpl2ldist, 1);
%         l2p(s,m).Complex.dist=mean(tmpl2pdist, 1);
%         p2p(s,m).Complex.dist=mean(tmpp2pdist, 1);
%         p2l(s,m).Complex.dist=mean(tmpp2ldist, 1);
%         clear tmpl2lacc tmpl2pacc tmpp2lacc tmpp2pacc tmpl2ldist tmpl2pdist tmpp2pdist tmpp2ldist
%         
%         l2l(s,m).Complex.thresh=quantile(l2l(s,m).Complex.dist, .95);
%         l2p(s,m).Complex.thresh=quantile(l2p(s,m).Complex.dist, .95);
%         p2p(s,m).Complex.thresh=quantile(p2p(s,m).Complex.dist, .95);
%         p2l(s,m).Complex.thresh=quantile(p2l(s,m).Complex.dist, .95);
%     end
%     waitbar(m/length(userOptions.maskNames));
%     
% end
% close(h)

save('Statistics/RBF_default_ROI_SVM.mat', 'l2l', 'l2p', 'p2p', 'p2l', ...
    'lComplex', 'pComplex', 'lanim', 'panim', 'lanimFam', 'panimFam', ...
    'linAnimFam', 'pinAnimFam', 'lfam', 'pfam', 'lverb', 'pverb', ...
    'larpc', 'parpc', 'lsyntax', 'psyntax', 'lrelcan', 'prelcan', ...
    'userOptions')
