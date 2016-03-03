%get similarity between stimuli
n=8;
s3={'the', 'gravy', 'touched', 'the', 'broccoli'};
s7={'the', 'trainer', 'stretched', 'the', 'dancer'};
s8={'the', 'candle', 'lit', 'the', 'match'};
s1={'the', 'bus', 'hit', 'the', 'truck'};
s6={'the', 'tree', 'crushed', 'the', 'car'};
s4={'the', 'goalie', 'kicked', 'the', 'referee'};
s2={'the', 'singer', 'kissed', 'the', 'guitarist'};
s5={'the', 'woman', 'consoled', 'the', 'man'};

smat={s1; s2; s3; s4; s5; s6; s7; s8};

% load('wiki/t300/wp_chain1_sample1.mat')
% load('wiki/wiki_vocab.mat')
%% load wp and vocab here %%
%% 
simMatrix=zeros(n);
for i = 1:length(smat)
    for j = 1:length(smat)
        simMatrix(i,j)=getSimilarity(smat{i}, smat{j}, wp, vocabulary, 'normal');
        simMatrix(j,i)=getSimilarity(smat{j}, smat{i}, wp, vocabulary, 'normal');
    end
end
simVec=reshape(simMatrix, 1, numel(simMatrix));
[~, ~, ranking]=unique(simMatrix);
simRank=reshape(ranking, n,n);

% save('simMatrix_v1.mat', 'simMatrix', 'simRank') %v1: using 300 topics
save('simMatrix_v2.mat', 'simMatrix', 'simRank') %2: using 1000 topics