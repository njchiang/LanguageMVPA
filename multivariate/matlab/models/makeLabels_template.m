function [ models ] = makeLabels_music()
%Making labels for SVM ROI Analysis
%   Detailed explanation goes here
    nConditions = 84;
    LMStimType=[ones(1, nConditions/2) 2*ones(1,nConditions/2)];
    models(1).label=LMStimType;
    models(1).name='LMStimType';
    LMSyntax=repmat([ones(1, nConditions/6) 2*ones(1,nConditions/6) 3*ones(1,nConditions/6)],1,2);
    models(2).label=LMSyntax;
    models(2).name='LMSyntax';
    LSyntax=LMSyntax(LMStimType==1);
    models(3).label=LSyntax;
    models(3).name='LSyntax';
    MSyntax=LMSyntax(LMStimType==2);
    models(4).label=MSyntax;
    models(4).name='MSyntax';
    LActPass=[ones(1,nConditions/6) 2*ones(1,nConditions/6) zeros(1,nConditions/6)];
    models(5).label=LActPass;
    models(5).name='LActPass';
    MActPass=LActPass;
    models(6).label=MActPass;
    models(6).name='MActPass';
    LvMSyntax=LSyntax;
    models(7).label=LvMSyntax;
    models(7).name='LvMSyntax';
    MvLSyntax=MSyntax;
    models(8).label=MvLSyntax;
    models(8).name='MvLSyntax';
end

