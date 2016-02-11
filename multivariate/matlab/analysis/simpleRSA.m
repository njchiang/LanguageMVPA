% take in RDM and extract between - within
% label=repmat([4 3 2 1], 1, 16)';
% label(33:64) = 0;
[diff, p] = simpleRSA(RDM, labels)

syntax=repmat([4 3 2 1], 1, 16)';
% actpass=repmat([2 2 1 1], 1, 16)';
% relcan=repmat([2 1 2 1], 1, 16)';
% verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,2)';
% stimType = [ones(1,32) 2*ones(1,32)];
% label = [1 2 1 2];
% label(3:4) = 0;


corrMat = 1-RDM;
withinMat = corrMat(label==1, label==1);
withinScore = mean(vectorizeRDM(withinMat));
nWithin = numel(vectorizeRDM(withinMat));

betweenMat = corrMat(label==1, label ~= 1 & label ~=0);
nBetween=numel(betweenMat);
betweenScore = mean(betweenMat(:));

[h, p] = ttest2(vectorizeRDM(withinMat)', reshape(betweenMat, numel(betweenMat),1));

diff = betweenScore - withinScore;

% maybe should do R to Z transform... then can just do a z test between the
% two?