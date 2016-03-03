function [topics p_topics excludedWords, p_z_w] = getTopicProbability(wp, vocabulary, inputWords)
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

%per Dawn's suggestion
p_w_z_all=bsxfun(@rdivide, wp,  sum(wp));
p_z_w_all=bsxfun(@rdivide, wp, sum(wp,2));
p_w=sum(prod(p_z_w_all(wordIdx,:)));

for t = 1:T
    p_z_w(t)=prod(p_z_w_all(wordIdx,t)*prior)/p_w;
end

topics=find(p_z_w >0);
[~, sortIdx]=sort(-p_z_w(p_z_w>0));
p_topics=p_z_w(p_z_w>0);
p_topics=p_topics(sortIdx);
topics=topics(sortIdx);

end


