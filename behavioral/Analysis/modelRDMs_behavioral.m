function Models = modelRDMs_behavioral()
genModel = ones(32,32);
syntaxLabels = repmat([1 2 3 4],1, 8);
Syntax=genModel;
for i = 1:max(syntaxLabels)
    Syntax(syntaxLabels == i, syntaxLabels ==i) = 0;
end
Models.Syntax = Syntax;
Semantics = genModel;
semanticLabels = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8];
for i = 1:max(semanticLabels)
    Semantics(semanticLabels == i, semanticLabels ==i) = 0;
end
Models.Semantics = Semantics;
end