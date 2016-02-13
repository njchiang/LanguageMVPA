try
    load('D:/CloudStation/Grad/Completed_Research/Topics/wiki/t300/wp_chain1_sample1.mat');
    load('D:/CloudStation/Grad/Completed_Research/Topics/wiki/wiki_vocab.mat');
catch
    load('~/CloudStation/Grad/Completed_Research/Topics/wiki/t300/wp_chain1_sample1.mat');
    load('~/CloudStation/Grad/Completed_Research/Topics/wiki/wiki_vocab.mat');
end
inputWords = {'broccoli', 'touched', 'gravy', 'candle', 'lit', 'match', ...
    'truck', 'hit', 'bus', 'tree', 'crushed', 'car', ...
    'singer', 'kissed', 'guitarist', 'dancer', 'stretched', 'trainer',  ...
    'goalie', 'kicked', 'referee', 'woman', 'consoled', 'man'};

[wordsGivenTopics, topicsGivenWords, excludedWords] = getTopicFeatureVectors(wp, vocabulary, inputWords);

goodTopics = (sum(wordsGivenTopics)~=0);
% filter out useless topics
singleWordRepres = wordsGivenTopics(:, goodTopics);
goodTopics = (sum(wordsGivenTopics>0)>1);
singleWordRepres = wordsGivenTopics(:,goodTopics);

singleWordDistn=bsxfun(@rdivide, singleWordRepres, sum(singleWordRepres,2)); % probability of word given topic
% filter out unique topics
% need to figure how to index properly, e.g. match index to correct
% regressor

for i = 1:8
    % which method to use?
    semanticsResSum(i,:) = sum(singleWordRepres(3*(i-1)+(1:3),:));
    semanticsResProd(i,:) = prod(singleWordRepres(3*(i-1)+(1:3),:));
end
semanticsResSumDistn=bsxfun(@rdivide, semanticsResSum, sum(semanticsResSum,2));
semanticsResProdDistn=bsxfun(@rdivide, semanticsResProd, sum(semanticsResProd,2));
csvwrite('D:\CloudStation\scratchpad\LanguageMVPA\topicRepresentation.csv', semanticsResSumDistn)
subplot(1,2,1)
imagesc(corr(semanticsResSumDistn')); colorbar;
subplot(1,2,2)
imagesc(corr(semanticsResProdDistn'));colorbar;

labels=kron([1 2], ones(1,4))';
libsvmtrain(labels, semanticsResSumDistn, '-s 0 -t 1 -v 4')

% get a cope for each topic

[COEFF, SCORE, LATENT] = pca(X)

X = semanticsResSumDistn;
% repeat each thing 4 times down
X = kron(X, ones(4,1));
Y = responsePatterns.left_IFG_operc.LMVPA005';
[B, FitInfo] = lassoglm(X, Y(1:32,1))
[B, FitInfo] = lasso(X, Y(1:32,1))

B = inv(X'*X)*(X'*Y(1:32,:))
%% I think i want to be working with words given topics

% goodTopics = (sum(topicsGivenWords)~=0);
% res = topicsGivenWords(:, goodTopics);