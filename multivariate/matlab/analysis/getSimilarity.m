% getSimilarity.m
% gets similarity between two word vectors in the topics model
% Griffiths, Steyvers, Tenenbaum 2007
function [condP, excludedWords1, excludedWords2] = getSimilarity(w1, w2, wp, vocabulary, analysisType)

% convert word vectors to index format
if strcmp(analysisType, 'log')
    disp('computing log likelihoods')
else
    disp('computing based on raw probabilities')
end
[wordIdx1, excludedWords1] = inputWordsToIdx(vocabulary, w1);
[wordIdx2, excludedWords2] = inputWordsToIdx(vocabulary, w2);
wp=full(wp);
p_w_z=bsxfun(@rdivide, wp,  sum(wp));

%move 0's to tiny number
p_w_z(p_w_z==0)=eps;
%from griffiths, steyvers, tenenbaum 2007
%conditional probability of w2|w1
% condP=1/(1+(sum((1-prod(p_w_z(wordIdx2,:))).*prod(p_w_z(wordIdx1,:)),2) ...
%     /sum(prod(p_w_z(wordIdx2,:)).*prod(p_w_z(wordIdx1,:)),2)));
% condP=sum(prod(p_w_z(wordIdx2,:),1).*prod(p_w_z(wordIdx1,:),1),2) ...
%     /sum(prod(p_w_z(wordIdx1,:),1),2);
%     /sum(prod(p_w_z(wordIdx2,:)).*prod(p_w_z(wordIdx1,:)),2)));

if strcmp(analysisType, 'log')
condP=sum(sum(log(p_w_z(wordIdx2,:)),1)+sum(log(p_w_z(wordIdx1,:)),1),2)...
    - sum(sum(log(p_w_z(wordIdx1,:)),1),2);
else
    condP=sum(prod(p_w_z(wordIdx2,:),1).*prod(p_w_z(wordIdx1,:),1),2) ...
    /sum(prod(p_w_z(wordIdx1,:),1),2);
end
end


