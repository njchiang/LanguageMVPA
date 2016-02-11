function Models = modelRDMs_avg_syntax()

syntax=repmat([4 3 2 1], 1, 2);
actpass=repmat([2 2 1 1], 1, 2);
relcan=repmat([2 1 2 1], 1, 2);

genModel=ones(8,8);
SyntaxMat=genModel;
SyntaxMat(syntax==1, syntax==1)=0;
SyntaxMat(syntax==2, syntax==2)=0;
SyntaxMat(syntax==3, syntax==3)=0;
SyntaxMat(syntax==4, syntax==4)=0;
SyntaxMat(logical(eye(size(genModel)))) = 0;
Models.SyntaxDetector=SyntaxMat;

ActPassMat=genModel;
ActPassMat(actpass==1, actpass==1)=0;
ActPassMat(actpass==2, actpass==2)=0;
ActPassMat(logical(eye(size(genModel))))=0;
RelCanMat=genModel;
RelCanMat(relcan==1, relcan==1)=0;
RelCanMat(relcan==2, relcan==2)=0;
RelCanMat(logical(eye(size(genModel))))=0;

Models.SyntaxComplex = ActPassMat+RelCanMat+SyntaxMat;

end