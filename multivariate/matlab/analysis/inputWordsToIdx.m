function [wordIdx, excludedWords] = inputWordsToIdx(vocabulary, inputWords)
excludedWords=[];
wordIdx=[];
nullCounter=0;
for i = 1:length(inputWords)
    if sum(strcmp(vocabulary, inputWords{i}))
        wordIdx=[wordIdx find(strcmp(vocabulary, inputWords{i}))];
    else
        nullCounter=nullCounter+1;
        excludedWords{nullCounter}=inputWords{i};
    end
end