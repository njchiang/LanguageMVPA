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
% only interested in topics in which more than one element is represented.
semanticsResSumDistn=bsxfun(@rdivide, semanticsResSum, sum(semanticsResSum,2));
semanticsResProdDistn=bsxfun(@rdivide, semanticsResProd, sum(semanticsResProd,2));

goodSumTopics = (sum(semanticsResSumDistn>0)>4);
mrmrFeat = [84 61 87 15 89 19 86 78 97 45 33 57 93 91 5 2 58 8]
semanticsResSumDistn2 = semanticsResSumDistn(:,mrmrFeat);
sum(semanticsResSumDistn2~=0,2)
animLabels=[1;1;1;1;2;2;2;2];
libsvmtrain(animLabels, semanticsResSumDistn2, '-s 0 -t 0 -v 2')
kmeans(semanticsResSumDistn2,2)
subplot(1,2,1)
imagesc(corr(semanticsResSumDistn2'))
subplot(1,2,2)
imagesc(corr(semanticsResSumDistn'))
csvwrite('D:\CloudStation\scratchpad\LanguageMVPA\topicRepresentation.csv', semanticsResSumDistn)


%% messing around
applyPCA=@(d, PCs)(d*PCs*inv(PCs'*PCs));
% doing algebra;
% testImages=score*pcs';
% testImages*pcs=score'*pcs*pcs
% testImages*pcs*inv(pcs'*pcs)=score*pcs'*pcs*inv(pcs'*pcs);
%% part 1
trainImages=rawImages(:,1:150)';
testImages=rawImages(:,151:end)';

% compute eigen faces based on training images
meanTrainImages=mean(trainImages);
centeredTrainImages=bsxfun(@minus, trainImages, meanTrainImages)';
%goal: eigenvectors of 1/M*centeredTrainImages*centeredTrainImages', but
%this is infeasible

[trainVecs, trainVals] = eig(centeredTrainImages'*centeredTrainImages); 
eigenFaces=centeredTrainImages*trainVecs;
dispEigenFaces=bsxfun(@plus, eigenFaces, meanTrainImages');

[pcs, score, latent]=pca(centeredTrainImages);

%center the test images
meanTestImages=mean(testImages);
centeredTestImages=bsxfun(@minus, testImages, meanTestImages);

%top 20 eigenfaces
figure;
for k = 1:150
    if k <= 20
    subplot(4,5,k)
    imagesc(reshape(eigenFaces(:,end-k+1), 256, 256)+.5)
    title(num2str(k))
    axis off
    end
    newScores=applyPCA(centeredTestImages, dispEigenFaces(:,1:k));
%     newImages=bsxfun(@plus,newScores*pcs(:,1:k)',meanTestImages);
    newImages=bsxfun(@plus,newScores*dispEigenFaces(:,1:k)',meanTestImages);
    % find the error between reconstructed (centered) image and centered
    % image
    totalReconError(k)=mean(((newImages(:)-testImages(:)).^2));
end

subplot(1,2,1)
imagesc(corr(semanticsResSumDistn')); colorbar;
subplot(1,2,2)
imagesc(corr(semanticsResProdDistn'));colorbar;

labels=kron([1 2], ones(1,4))';
labels = repmat([4;3;2;1], 8,1);
labels = repmat([3;2;2;1], 8,1);
labels = kron([8;7;6;5;4;3;2;1], ones(4,1));
labels = kron([1;2;3;4;5;6;7;8], ones(4,1));

libsvmtrain(labels, Y, '-s 1 -t 0 -v 4')
libsvmtrain(labels, Y, '-s 4 -t 0 -v 4')

% get a cope for each topic

[COEFF, SCORE, LATENT] = pca(X)

X = semanticsResSumDistn;
% repeat each thing 4 times down
% X = semanticsResSumDistn2;
X = kron(X, ones(4,1));
Y = responsePatterns.left_IFG_operc.LMVPA005(:,1:32)';
[B,STATS] = lassoglm(X,Y)
[B, FitInfo] = lassoglm(X, Y(1:32,1))
[B, FitInfo] = lasso(X, Y(1:32,1), 'CV', 4,'Alpha',.5);
lassoPlot(B,FitInfo,'plottype','CV')
B = inv(X'*X)*(X'*Y(1:32,:))

% playing with mitchell data
idx1=[];idx2=[];
for i = 1:length(info)
if strcmp(info(i).word, 'saw'); idx1 = [idx1 i];  end
if strcmp(info(i).word, 'cow'); idx2 = [idx2 i]; end
%% seems to be able to do different categories...
% if strcmp(info(i).cond, 'manmade'); idx1 = [idx1 i];  end
% if strcmp(info(i).cond, 'vehicle'); idx2 = [idx2 i]; end
end
labels1 = ones(size(idx1));
labels2 = 2*ones(size(idx2));
idx = [idx1 idx2]';
labels = [labels1 labels2]';
libsvmtrain(labels, mitchellData(idx,:), '-v 2 -s 0 -t 0')
imagesc(1-corr(mitchellData(idx,:)'))
%% I think i want to be working with words given topics

% goodTopics = (sum(topicsGivenWords)~=0);
% res = topicsGivenWords(:, goodTopics);