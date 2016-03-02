%  modelRDMs is a user-editable function which specifies the models which
%  brain-region RDMs should be compared to, and which specifies which kinds of
%  analysis should be performed.
%
%  Models should be stored in the "Models" struct as a single field labeled
%  with the model's name (use underscores in stead of spaces).
%
%  Cai Wingfield 11-2009

function Models = modelRDMs_category_ROI()
anim=repmat([1*ones(1,16) 2*ones(1,16)], 1, 2);
% fam=repmat([1*ones(1,4) 2*ones(1,8) 1*ones(1,4) 3*ones(1,8) 4*ones(1,8) ], 1, 2);
verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,2);
syntax=repmat([4 3 2 1], 1, 16);
actpass=repmat([2 2 1 1], 1, 16);
relcan=repmat([2 1 2 1], 1, 16);
complexlab=repmat([2 2 2 1], 1, 16);


genModel=ones(64,64);
SyntaxMat=genModel;
SyntaxMat(syntax==1, syntax==1)=0;
SyntaxMat(syntax==2, syntax==2)=0;
SyntaxMat(syntax==3, syntax==3)=0;
SyntaxMat(syntax==4, syntax==4)=0;
SyntaxMat(logical(eye(size(genModel)))) = 0;

[Models.LSyntaxDetector, Models.PSyntaxDetector, Models.CSyntaxDetector] = ...
    splitModels(SyntaxMat);

ActPassMat=genModel;
ActPassMat(actpass==1, actpass==1)=0;
ActPassMat(actpass==2, actpass==2)=0;
ActPassMat(logical(eye(size(genModel))))=0;

[Models.LActPassDetector, Models.PActPassDetector, Models.CActPassDetector] = ...
    splitModels(ActPassMat);

RelCanMat=genModel;
RelCanMat(relcan==1, relcan==1)=0;
RelCanMat(relcan==2, relcan==2)=0;
RelCanMat(logical(eye(size(genModel))))=0;

[Models.LRelCanDetector, Models.PRelCanDetector, Models.CRelCanDetector] = ...
    splitModels(RelCanMat);

SyntaxComplexityMat = ActPassMat + RelCanMat;
[Models.LComplexity, Models.PComplexity, Models.CComplexity] = ...
    splitModels(SyntaxComplexityMat);

VerbMat=genModel;
VerbMat(verb==1, verb==1)=0;
VerbMat(verb==2, verb==2)=0;
VerbMat(verb==3, verb==3)=0;
VerbMat(verb==4, verb==4)=0;
VerbMat(verb==5, verb==5)=0;
VerbMat(verb==6, verb==6)=0;
VerbMat(verb==7, verb==7)=0;
VerbMat(verb==8, verb==8)=0;
VerbMat(logical(eye(size(genModel))))=0;

[Models.LVerbDetector, Models.PVerbDetector, Models.CVerbDetector] = ...
    splitModels(RelCanMat);

AnimMat=genModel;
AnimMat(anim==1, anim==1)=0;
AnimMat(anim==2, anim==2)=0;

[Models.LAnimDetector, Models.PAnimDetector, Models.CAnimDetector] = ...
    splitModels(AnimMat);

VerbHierarchy = AnimMat + VerbMat;
[Models.LVerbModel, Models.PVerbModel, Models.CVerbModel] = ...
    splitModels(VerbHierarchy);

end

function [L, P, C] = splitModels(orig)
L = orig;
L(33:64,:) = nan;
L(:, 33:64)=nan;
L(logical(eye(size(orig))))=0;
P = orig;
P(1:32, :) = nan;
P(:, 1:32) = nan;
P(logical(eye(size(orig))))=0;
C = orig;
C(1:32, 1:32) = nan;
C(33:64, 33:64) = nan;
C(logical(eye(size(orig))))=0;
end