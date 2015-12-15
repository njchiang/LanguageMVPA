%  modelRDMs is a user-editable function which specifies the models which
%  brain-region RDMs should be compared to, and which specifies which kinds of
%  analysis should be performed.
%
%  Models should be stored in the "Models" struct as a single field labeled
%  with the model's name (use underscores in stead of spaces).
%
%  Cai Wingfield 11-2009

function Models = modelRDMs_LSA_L()
anim=repmat([1*ones(1,16) 2*ones(1,16)], 1, 2);
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

Models.LSyntax=SyntaxMat;
Models.LSyntax(33:64, 33:64)=nan;
Models.LSyntax(logical(eye(size(genModel))))=0;

ActPassMat=genModel;
ActPassMat(actpass==1, actpass==1)=0;
ActPassMat(actpass==2, actpass==2)=0;
ActPassMat(logical(eye(size(genModel))))=0;

Models.LActPass=ActPassMat;
Models.LActPass(33:64, 33:64)=nan;
Models.LActPass(logical(eye(size(genModel))))=0;

RelCanMat=genModel;
RelCanMat(relcan==1, relcan==1)=0;
RelCanMat(relcan==2, relcan==2)=0;
RelCanMat(logical(eye(size(genModel))))=0;

Models.LRelCan=RelCanMat;
Models.LRelCan(33:64, 33:64)=nan;
Models.LRelCan(logical(eye(size(genModel))))=0;

StructureMat=genModel;
% StructureMat(complexlab==1, complexlab==1)=0;
StructureMat(complexlab==2, complexlab==2)=0;

Models.LStructure=StructureMat;
Models.LStructure(33:64, 33:64)=nan;
Models.LStructure(logical(eye(size(genModel))))=0;

Models.LCSyntax=Models.LSyntax+Models.LRelCan+Models.LActPass; %+Models.LStructure;

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

Models.LVerb=VerbMat;
Models.LVerb(33:64, 33:64)=nan;
Models.LVerb(logical(eye(size(genModel))))=0;

AnimMat=genModel;
AnimMat(anim==1, anim==1)=0;
AnimMat(anim==2, anim==2)=0;

Models.LAnim=AnimMat;
Models.LAnim(33:64, 33:64)=nan;
Models.LAnim(logical(eye(size(genModel))))=0;

Models.LSemantics=Models.LVerb+Models.LAnim;