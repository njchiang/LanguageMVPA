load('tstat_responsePatterns.mat')

syntax=repmat([4 3 2 1], 1, 8)';
verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,1)';

labels = syntax;
prior = [10 5 3 1];
runs = verb;

labels = verb;
prior = [1 1 1 1 1 1 1 1];
runs = syntax;

masks = fieldnames(responsePatterns);
subs = fieldnames(responsePatterns.(masks{1}));
gnbAccs = zeros(length(subs), length(masks)-1);
mapAccs = zeros(length(subs), length(masks)-1);
svmAccs = zeros(length(subs), length(masks)-1);
for m = 1:length(masks)-1
    for s = 1:length(subs)
        pat = responsePatterns.(masks{m}).(subs{s});
        pat = pat(:,1:32)';
        pat = zscore(pat);
        prior = prior ./sum(prior); % make sure sum to 1
        nFolds = length(unique(runs));
        classLab = zeros(nFolds, length(unique(trainLabs)));
        MAP = zeros(nFolds, length(unique(trainLabs)));
        gnbAcc = zeros(nFolds,1);
        mapAcc = zeros(nFolds,1);
        svmAcc = zeros(nFolds,1);
        for f = 1:nFolds
            trainSet = pat(runs(1:32)~=f,:);
            testSet = pat(runs(1:32)==f,:);
            
            trainLabs = labels(runs~=f,:);
            testLabs = labels(runs==f,:);
            imDatabase = zeros(length(unique(trainLabs)), size(pat,2));
            for i = 1:length(unique(trainLabs))
                imDatabase(i,:) = mean(trainSet(trainLabs==i,:));
            end
            
            classOut = corr(imDatabase', testSet'); % each column is a set of decision values
            [~, classLab(f,:)] = max(classOut);
            normedClassOut = bsxfun(@times, classOut, 1./sum(classOut));
            posteriorOut = bsxfun(@times, normedClassOut, prior);
            [~, MAP(f,:)] = max(posteriorOut);
            gnbAcc(f) = sum(classLab(f,:)' == testLabs) / length(testLabs);
            mapAcc(f) = sum(MAP(f,:)' == testLabs) / length(testLabs);
            svmStruct=libsvmtrain(trainLabs, trainSet, '-s 0 -t 1');
            [predicted_label, accuracy, ~] = libsvmpredict(testLabs, testSet, svmStruct);
            svmAcc(f) = accuracy(1);
        end
        
        gnbAccs(s,m) = mean(gnbAcc);
        mapAccs(s,m) = mean(mapAcc);
        svmAccs(s,m) = mean(svmAcc)/100;
    end
end