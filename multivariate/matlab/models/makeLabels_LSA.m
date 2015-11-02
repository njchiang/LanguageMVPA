function [ models ] = makeLabels_LSA()
%Making labels for SVM ROI Analysis
%   Detailed explanation goes here
nConditions = 64;
StimType=[ones(1,nConditions/2) 2*ones(1, nConditions/2)]';
%1: inanimate, 2: animate
anim=repmat([1*ones(1,4) 2*ones(1,4) 1*ones(1,4) 2*ones(1,8) 1*ones(1,8) 2*ones(1,4)], 1, 2)';
fam=repmat([1*ones(1,4) 2*ones(1,8) 1*ones(1,4) 3*ones(1,8) 4*ones(1,8) ], 1, 2)';
verb=repmat([ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)],1,2)';
syntax=repmat([4 3 2 1], 1, 16)';
actpass=repmat([2 2 1 1], 1, 16)';
relcan=repmat([2 1 2 1], 1, 16)';
complex=repmat([2 2 2 1], 1, 16)';
ARPC=repmat([0 2 1 0], 1, 16)';

models(1).label=StimType;
models(1).name='StimType';
models(1).runs=syntax;
% 8 way verb classification
LVerb=verb;
LVerb(StimType==2)=0;
models(2).label=LVerb;
models(2).name='LVerb';
models(2).runs=syntax;
% 2 way anim classification
LAnim=anim;
LAnim(StimType==2)=0;
models(3).label=LAnim;
models(3).name='LAnim';
models(3).runs=syntax;

% Animate: family
LAmFam=LVerb;
LAmFam(LAnim==1)=0;
models(4).label=LAmFam;
models(4).name='LAmFam';
models(4).runs=syntax;

LInAmFam=LVerb;
LInAmFam(LAnim==2)=0;
models(5).label=LInAmFam;
models(5).name='LInAmFam';
models(5).runs=syntax;


% 8 way verb classification
PVerb=verb;
PVerb(StimType==1)=0;
models(6).label=PVerb;
models(6).name='PVerb';
models(6).runs=syntax;

% 2 way anim classification
PAnim=anim;
PAnim(StimType==1)=0;
models(7).label=PAnim;
models(7).name='PAnim';
models(7).runs=syntax;

% Animate: family
PAmFam=PVerb;
PAmFam(PAnim==1)=0;
models(8).label=PAmFam;
models(8).name='PAmFam';
models(8).runs=syntax;

PInAmFam=PVerb;
PInAmFam(PAnim==2)=0;
models(9).label=PInAmFam;
models(9).name='PInAmFam';
models(9).runs=syntax;

% Syntax: 4-way
LSyntax=syntax;
LSyntax(StimType==2)=0;
models(10).label=LSyntax;
models(10).name='LSyntax';
models(10).runs=verb;

LActPass=actpass;
LActPass(StimType==2)=0;
models(11).label=LActPass;
models(11).name='LActPass';
models(11).runs=verb;

% Relative vs Canonical
LRelCan=relcan;
LRelCan(StimType==2)=0;
models(12).label=LRelCan;
models(12).name='LRelCan';
models(12).runs=verb;


% AR vs PC
LARPC=ARPC;
LARPC(StimType==2)=0;
models(13).label=LARPC;
models(13).name='LARPC';
models(13).runs=verb;

PSyntax=syntax;
PSyntax(StimType==1)=0;
models(14).label=PSyntax;
models(14).name='PSyntax';
models(14).runs=verb;

PActPass=actpass;
PActPass(StimType==1)=0;
models(15).label=PActPass;
models(15).name='PActPass';
models(15).runs=verb;

% Relative vs Canonical
PRelCan=relcan;
PRelCan(StimType==1)=0;
models(16).label=PRelCan;
models(16).name='PRelCan';
models(16).runs=verb;

% AR vs PC
PARPC=ARPC;
PARPC(StimType==1)=0;
models(17).label=PARPC;
models(17).name='PARPC';
models(17).runs=verb;
% % 4 way family classification
% LFam=fam;
% LFam(StimType==2)=0;
% models(3).label=LFam;
% models(3).name='LFam';
% models(3).runs=syntax;
% 
% % 2 way anim classification
% LAnim=anim;
% LAnim(StimType==2)=0;
% models(4).label=LAnim;
% models(4).name='LAnim';
% models(4).runs=syntax;
% 
% % Animate: family
% LAmFam=LFam;
% LAmFam(LAnim==1)=0;
% models(5).label=LAmFam;
% models(5).name='LAmFam';
% models(5).runs=syntax;
% 
% LInAmFam=LFam;
% LInAmFam(LAnim==2)=0;
% models(6).label=LInAmFam;
% models(6).name='LInAmFam';
% models(6).runs=syntax;

% % 8 way verb classification
% PVerb=verb;
% PVerb(StimType==1)=0;
% models(7).label=PVerb;
% models(7).name='PVerb';
% models(7).runs=syntax;
% % 
% % 4 way family classification
% PFam=fam;
% PFam(StimType==1)=0;
% models(8).label=PFam;
% models(8).name='PFam';
% models(8).runs=syntax;
% 
% % 2 way anim classification
% PAnim=anim;
% PAnim(StimType==1)=0;
% models(9).label=PAnim;
% models(9).name='PAnim';
% models(9).runs=syntax;
% 
% % Animate: family
% PAmFam=PFam;
% PAmFam(PAnim==1)=0;
% models(10).label=PAmFam;
% models(10).name='PAmFam';
% models(10).runs=syntax;
% 
% PInAmFam=PFam;
% PInAmFam(PAnim==2)=0;
% models(11).label=PInAmFam;
% models(11).name='PInAmFam';
% models(11).runs=syntax;
% 
% % Syntax: 4-way
% LSyntax=syntax;
% LSyntax(StimType==2)=0;
% models(12).label=LSyntax;
% models(12).name='LSyntax';
% models(12).runs=verb;
% 
% LActPass=actpass;
% LActPass(StimType==2)=0;
% models(13).label=LActPass;
% models(13).name='LActPass';
% models(13).runs=verb;
% 
% % Relative vs Canonical
% LRelCan=relcan;
% LRelCan(StimType==2)=0;
% models(14).label=LRelCan;
% models(14).name='LRelCan';
% models(14).runs=verb;
% 
% % Complex vs Simple
% LComplex=complex;
% LComplex(StimType==2)=0;
% models(15).label=LComplex;
% models(15).name='LComplex';
% models(15).runs=verb;
% 
% % AR vs PC
% LARPC=ARPC;
% LARPC(StimType==2)=0;
% models(16).label=LARPC;
% models(16).name='LARPC';
% models(16).runs=verb;
% 
% PSyntax=syntax;
% PSyntax(StimType==1)=0;
% models(17).label=PSyntax;
% models(17).name='PSyntax';
% models(17).runs=verb;
% 
% PActPass=actpass;
% PActPass(StimType==1)=0;
% models(18).label=PActPass;
% models(18).name='PActPass';
% models(18).runs=verb;
% 
% % Relative vs Canonical
% PRelCan=relcan;
% PRelCan(StimType==1)=0;
% models(19).label=PRelCan;
% models(19).name='PRelCan';
% models(19).runs=verb;
% 
% % Complex vs Simple
% PComplex=complex;
% PComplex(StimType==1)=0;
% models(20).label=PComplex;
% models(20).name='PComplex';
% models(20).runs=verb;
% 
% % AR vs PC
% PARPC=ARPC;
% PARPC(StimType==1)=0;
% models(21).label=PARPC;
% models(21).name='PARPC';
% models(21).runs=verb;

end

