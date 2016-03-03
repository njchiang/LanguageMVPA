function [wordsGivenTopics, topicsGivenWords, excludedWords] = getTopicFeatureVectors(wp, vocabulary, inputWords)
% rewrite this to return the whole vector as well as a mask
%wp is the output of topics model
% this calculates the probability for each of these words...
% output: topic idx, from greatest to lowest
% model hyperparameters

T=size(wp, 2); %T topics
wp=full(wp);
prior=1/T;
BETA = 200 / length(vocabulary);%0.01;
ALPHA = 50 / T;%0.1;


[wordIdx, excludedWords] = inputWordsToIdx(vocabulary, inputWords);

%per Dawn's suggestion... but why convert to posterior?
%wp: rows: words cols: topics
% probability of word given topic
p_w_z_all=bsxfun(@rdivide, wp,  sum(wp)); % sum for same topic over all words
% probability of topic given word
p_z_w_all=bsxfun(@rdivide, wp, sum(wp,2)); % sum for same word over all topics

% not sure which one to use for now... but they're highly correlated. add
% these as covariates after sparisfying them?

topicsGivenWords=p_z_w_all(wordIdx,:);

wordsGivenTopics = p_w_z_all(wordIdx,:);

% I just need to take the 300 vectors and use them as a feature...


end

