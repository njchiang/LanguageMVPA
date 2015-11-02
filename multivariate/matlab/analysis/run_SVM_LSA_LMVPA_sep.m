function [l2l, l2p, p2p, p2l, l_opt_params, p_opt_params] = run_SVM_LSA_LMVPA(labels, permutedLabels, currPat, runs, svmKernel)
% separate the correct voxels
sentPats=currPat(:,1:32)';
lPats=sentPats(labels>0,:);
picPats=currPat(:,33:64)';
pPats=picPats(labels>0,:);
numTrials=size(lPats,1);
labels=permutedLabels(:,1);
condLabels=unique(labels);
numPerCond=length(labels)/length(unique(labels));
for i = 1:10
    %need to take out half of each condition...0
    trainIND=[];
    testIND=[];
    for j=1:length(condLabels)
        b=condLabels(j);
        randIND=randperm(numPerCond);
        condIdx=find(labels==j);  
        trainIND=[trainIND; condIdx(randIND(1:numPerCond/2))];
        testIND=[testIND; condIdx(randIND(numPerCond/2+1:end))];
    end
    
%     trainIND=randIND(1:(numTrials/2)); % randomly choose half of the trials for optimization
%     testIND=randIND((numTrials/2+1):end);
    % optimize on random half of data
    cv_l_opt_params=optimizeSVM(permutedLabels(trainIND,1), lPats(trainIND,:), svmKernel);
    if svmKernel==2
        opts=['-s 0 -t ' num2str(svmKernel) ' -c ' num2str(cv_l_opt_params.best_C) ' -g ' num2str(cv_l_opt_params.best_gamma)];
    else
        opts=['-s 1 -t ' num2str(svmKernel) ' -n ' num2str(cv_l_opt_params.best_C) ];
    end
%     lCVStruct=libsvmtrain(permutedLabels(trainIND,1), lPats(trainIND,:), opts);
%     [l2l.acc(i), l2l.dist(i,:), l2l.thresh(i)]=permutationSVM(lCVStruct, lPats(testIND,:), permutedLabels(testIND,:), runs(testIND));
    
    %     [l2l.acc(i), l2l.dist(i,:), l2l.thresh(i)]=permutationSVM(lPats, permutedLabels, runs, opts);
    [l2l.acc(i), l2l.dist(i,:), l2l.thresh(i)] = CVpermutationSVM(lPats(testIND,:), permutedLabels(testIND,:), runs(testIND), opts);
    [l2l.acc(i), l2l.dist(i,:), l2l.thresh(i)] = CVpermutationSVM(lPats, permutedLabels, runs, opts);
    
    cv_p_opt_params=optimizeSVM(permutedLabels(trainIND,1), pPats(trainIND,:), svmKernel);
    if svmKernel==2
        opts=['-s 0 -t ' num2str(svmKernel) ' -c ' num2str(cv_p_opt_params.best_C) ' -g ' num2str(cv_p_opt_params.best_gamma)];
    else
        opts=['-s 1 -t ' num2str(svmKernel) ' -n ' num2str(cv_p_opt_params.best_C) ];
    end
%     pCVStruct=libsvmtrain(permutedLabels(trainIND,1), pPats(trainIND,:), opts);
%     [p2p.acc(i), p2p.dist(i,:), p2p.thresh(i)]=permutationSVM(pCVStruct, pPats(testIND,:), permutedLabels(testIND,:), runs(testIND));
    [p2p.acc(i), p2p.dist(i,:), p2p.thresh(i)] = CVpermutationSVM(pPats(testIND,:), permutedLabels(testIND,:), runs(testIND), opts);

end

l2l.acc=mean(l2l.acc);
l2l.dist=mean(l2l.dist);
l2l.thresh=quantile(l2l.dist,.95);
p2p.acc=mean(p2p.acc);
p2p.dist=mean(p2p.dist);
p2p.thresh=mean(p2p.dist, .95);

p_opt_params=optimizeSVM(permutedLabels(:,1), pPats, svmKernel);
if svmKernel==2
    opts=['-s 0 -t ' num2str(svmKernel) ' -c ' num2str(p_opt_params.best_C) ' -g ' num2str(p_opt_params.best_gamma)];
else
    opts=['-s 1 -t ' num2str(svmKernel) ' -n ' num2str(p_opt_params.best_C) ];
end
pStruct=libsvmtrain(permutedLabels(:,1), pPats, opts);
[p2l.acc, p2l.dist, p2l.thresh] = permutationSVM(pStruct, lPats, permutedLabels, runs);

l_opt_params=optimizeSVM(permutedLabels(:,1), lPats, svmKernel);
if svmKernel==2
    opts=['-s 0 -t ' num2str(svmKernel) ' -c ' num2str(l_opt_params.best_C) ' -g ' num2str(l_opt_params.best_gamma)];
else
    opts=['-s 1 -t ' num2str(svmKernel) ' -n ' num2str(l_opt_params.best_C) ];
end
lStruct=libsvmtrain(permutedLabels(:,1), lPats, opts);
[l2p.acc, l2p.dist, l2p.thresh]=permutationSVM(lStruct, pPats, permutedLabels, runs);
