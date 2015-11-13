function [l2l, l2p, p2p, p2l, l_opt_params, p_opt_params] = run_SVM_LSA_LMVPA(labels, permutedLabels, currPat, runs, svmKernel)
% define function
applyZ=@(betas, mu, sigma)(bsxfun(@times, bsxfun(@minus, betas, mu), 1./sigma));
% separate the correct voxels
sentPats=currPat(1:32,:);
% sentPats=zscore(currPat(1:32,:));
lPats=sentPats(labels>0,:);
picPats=currPat(33:64,:);
% picPats=zscore(currPat(33:64,:));
pPats=picPats(labels>0,:);
numTrials=size(lPats,1);
labels=permutedLabels(:,1);
condLabels=unique(labels);
numPerCond=length(labels)/length(unique(labels));
if svmKernel==2
    opts=['-q -s 0 -t ' num2str(svmKernel)];
else
    opts=['-q -s 1 -t ' num2str(svmKernel)];
end
[l2l.acc, l2l.dist, l2l.thresh]=CVpermutationSVM(lPats, permutedLabels, runs, 'opts', opts);
%     cv_p_opt_params=optimizeSVM(permutedLabels(trainIND,1), pPats(trainIND,:), svmKernel);
if svmKernel==2
    opts=['-q -s 0 -t ' num2str(svmKernel)];
else
    opts=['-q -s 1 -t ' num2str(svmKernel)];
end
[p2p.acc, p2p.dist, p2p.thresh] = CVpermutationSVM(pPats, permutedLabels, runs, 'opts', opts);

[z_pPats, p_mu, p_sigma] = zscore(zscore(pPats,0,2));
plPats=applyZ(zscore(lPats,0,2), p_mu,1./p_sigma);
p_opt_params=optimizeSVM(permutedLabels(:,1), z_pPats, svmKernel);
if svmKernel==2
    opts=['-q -s 0 -t ' num2str(svmKernel) ' -c ' num2str(p_opt_params.best_C) ' -g ' num2str(p_opt_params.best_gamma)];
else
    opts=['-q -s 1 -t ' num2str(svmKernel) ' -n ' num2str(p_opt_params.best_C) ];
end
pStruct=libsvmtrain(permutedLabels(:,1), z_pPats, opts);
[p2l.acc, p2l.dist, p2l.thresh] = permutationSVM(pStruct, plPats, permutedLabels, runs);

[z_lPats, l_mu, l_sigma] = zscore(zscore(lPats,0,2));
lpPats=applyZ(zscore(pPats,0,2), l_mu, l_sigma);
l_opt_params=optimizeSVM(permutedLabels(:,1), z_lPats, svmKernel);
if svmKernel==2
    opts=['-q -s 0 -t ' num2str(svmKernel) ' -c ' num2str(l_opt_params.best_C) ' -g ' num2str(l_opt_params.best_gamma)];
else
    opts=['-q -s 1 -t ' num2str(svmKernel) ' -n ' num2str(l_opt_params.best_C) ];
end
lStruct=libsvmtrain(permutedLabels(:,1), z_lPats,  opts);
[l2p.acc, l2p.dist, l2p.thresh]=permutationSVM(lStruct, lpPats, permutedLabels, runs);
