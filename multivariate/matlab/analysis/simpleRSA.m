% take in RDM and extract between - within
% label=repmat([4 3 2 1], 1, 16)';
% label(33:64) = 0;
% syntax=repmat([4 3 2 1], 1, 16)';
% lsyntax = syntax;
% lsyntax(33:64) = 0;
% % actpass=repmat([2 2 1 1], 1, 16)';
% % relcan=repmat([2 1 2 1], 1, 16)';
% % verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,2)';
% % stimType = [ones(1,32) 2*ones(1,32)];
% % label = [1 2 1 2];
% % label(3:4) = 0;
% 
% for i = 1:16
%     thisRDM = RDMs(4,i).RDM;
%     [b(i),w(i)] = simpleRSA(thisRDM, lsyntax);
% end


function [betweenScore, withinScore, p] = simpleRSA(RDM, label)

% syntax=repmat([4 3 2 1], 1, 16)';
% actpass=repmat([2 2 1 1], 1, 16)';
% relcan=repmat([2 1 2 1], 1, 16)';
% verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,2)';
% stimType = [ones(1,32) 2*ones(1,32)];
% label = [1 2 1 2];
% label(3:4) = 0;
corrMat = 1-RDM;

withinScores = zeros(numel(vectorizeRDM(corrMat(label==1, label==1))),max(label));
betweenScores = zeros(numel(corrMat(label==1, label~=1 & label ~= 0)),max(label));
for i = 1:max(label)
    withinMat = corrMat(label==i, label==i);
    withinScores(:,i) =(fisherTransform(vectorizeRDM(withinMat)));
    nWithin = numel(vectorizeRDM(withinMat));
    
    betweenMat = corrMat(label==i, label ~= i & label ~=0);
    nBetween=numel(betweenMat);
    betweenScores(:,i) = reshape(fisherTransform(betweenMat(:)),nBetween,1);
end


[h, p] = ttest2(reshape(withinScores, numel(withinScores),1), reshape(betweenScores, numel(betweenScores),1));
betweenScore = mean(betweenScores(:));
withinScore = mean(withinScores(:));
diff = betweenScore - withinScore;
% maybe should do R to Z transform... then can just do a z test between the
% two?