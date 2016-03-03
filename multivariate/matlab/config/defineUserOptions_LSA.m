function userOptions = defineUserOptions_LSA()
%
%  projectOptions is a nullary function which initialises a struct
%  containing the preferences and details for a particular project.
%  It should be edited to taste before a project is run, and a new
%  one created for each substantially different project (though the
%  options struct will be saved each time the project is run under
%  a new name, so all will not be lost if you don't do this).
%
%  For a guide to how to fill out the fields in this file, consult
%  the documentation folder (particularly the userOptionguide.m)
%
%  Cai Wingfield 11-2009
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council


%% Project details

% This name identifies a collection of files which all belong to the same run of a project.
userOptions.analysisName = 'Compute';

% This is the root directory of the project.
% userOptions.rootPath = '/space/raid5/data/monti/Analysis/Music/mvpa_data';
userOptions.rootPath=pwd;
% The path leading to where the scans are stored (not including subject-specific identifiers).
% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
% "[[betaIdentifier]]" should be used as a placeholder to denote an output of betaCorrespondence.m if SPM is not being used; or an arbitrary filename if SPM is being used.
userOptions.betaPath = [userOptions.rootPath '/data/betas/LSA/[[subjectName]]_[[betaIdentifier]]'];% e.g. /imaging/mb01/lexpro/multivariate/ffx_simple/[[subjectName]]/[[betaIdentifier]]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEATUERS OF INTEREST SELECTION OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% %% %% %%
%% fMRI  %% Use these next three options if you're working in fMRI native space:
%% %% %% %% %%

% The path to a stereotypical mask data file is stored (not including subject-specific identifiers).
% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
% "[[maskName]]" should be used as a placeholder to denote an entry in userOptions.maskNames
userOptions.maskPath = [userOptions.rootPath '/data/masks/3mm_[[maskName]].img'];%'/imaging/mb01/lexpro/multivariate/ffx_simple/[[subjectName]]/[[maskName]].img';

% The list of mask filenames (minus .hdr extension) to be used.
userOptions.maskNames = { ...
    'left_IFG_operc',  ...
    'left_IFG_triang', ...
    'left_IFG_orbit', ...
    'left_STG_post',  ...
    'left_STG_ant', ...
    'left_SPL', ...
    'left_angular',  ...
    'left_MTG_post',...
    'left_caudate', ...
    'right_IFG_operc',  ...
    'right_IFG_triang', ...
    'right_IFG_orbit', ...
    'right_STG_post',  ...
    'right_STG_ant', ...
    'right_SPL', ...
    'right_angular',  ...
    'right_MTG_post',...
    'right_caudate', ...
    'grayMatter' ...
    };
%wholehead

% 		userOptions.maskNames = { 'wholehead' };
% userOptions.maskNames = { ...
% 'left_angular',     'left_MTG_post', 'right_BA47pos'...
% 'left_BA44_d',     'left_SPL',      'right_caudate',...
% 'left_BA44',       'left_STG_ant',  'right_IFG_operc',...
% 'left_BA44_v',      'left_STG_post',  'right_IFG_triang', ...
% 'left_BA45ant',     'right_angular',  'right_MTG_post', ...
% 'left_BA47ant',     'right_BA44_d',   'right_SPL', ...
% 'left_BA47',        'right_BA44',     'right_STG_ant', ...
% 'left_BA47pos',     'right_BA44_v',   'right_STG_post', ...
% 'left_caudate',     'right_BA45ant',   ...
% 'left_IFG_operc',   'right_BA47ant',  ...
% 'left_IFG_triang',  'right_BA47', 'wholehead'};

%%%%%%%%%%%%%%%%%%%%%%%%%
%% SEARCHLIGHT OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% %% %% %%
%% fMRI  %% Use these next three options if you're working in fMRI native space:
%% %% %% %% %%

% What is the path to the anatomical (structural) fMRI scans for each subject?
% "[[subjectName]]" should be used to denote an entry in userOptions.subjectNames
%NOT BEING USED
% userOptions.structuralsPath = [userOptions.rootPath '/data/MNI152_T1_3mm'];% e.g. /imaging/mb01/lexpro/[[subjectName]]/structurals/

% What are the dimensions (in mm) of the voxels in the scans?
userOptions.voxelSize = [3 3 3.75];

% What radius of searchlight should be used (mm)?
userOptions.searchlightRadius = 12;

%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL SETUP %%
%%%%%%%%%%%%%%%%%%%%%%%%

% The list of subjects to be included in the study.
userOptions.subjectNames = { ...
    'LMVPA001', 'LMVPA002', 'LMVPA003', 'LMVPA005', ...
    'LMVPA006', 'LMVPA007', 'LMVPA008', 'LMVPA009', ...
    'LMVPA010', 'LMVPA011', 'LMVPA013', 'LMVPA014', ...
    'LMVPA015', 'LMVPA016', 'LMVPA017', 'LMVPA018', ...
    'LMVPA019'};% eg CBUXXXXX LMVPA007', 'LMVPA001', 'LMVPA002', 
% The default colour label for RDMs corresponding to RoI masks (as opposed to models).
userOptions.RoIColor = [0 0 1];
userOptions.ModelColor = [0 1 0];

% Should information about the experimental design be automatically acquired from SPM metadata?
% If this option is set to true, the entries in userOptions.conditionLabels MUST correspond to the names of the conditions as specified in SPM.
userOptions.getSPMData = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS PREFERENCES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First-order analysis

% Text lables which may be attached to the conditions for MDS plots.
userOptions.conditionLabels={ ...
    'Touch_4', ...
    'Touch_3', ...
    'Touch_2', ...
    'Touch_1', ...
    'Stretch_4', ...
    'Stretch_3', ...
    'Stretch_2', ...
    'Stretch_1', ...
    'Light_4', ...
    'Light_3', ...
    'Light_2', ...
    'Light_1', ...
    'Kiss_4', ...
    'Kiss_3', ...
    'Kiss_2', ...
    'Kiss_1', ...
    'Kick_4', ...
    'Kick_3', ...
    'Kick_2', ...
    'Kick_1', ...
    'Hit_4', ...
    'Hit_3', ...
    'Hit_2', ...
    'Hit_1', ...
    'Crush_4', ...
    'Crush_3', ...
    'Crush_2', ...
    'Crush_1', ...
    'Console_4', ...
    'Console_3', ...
    'Console_2', ...
    'Console_1', ...
    'Touch_4', ...
    'Touch_3', ...
    'Touch_2', ...
    'Touch_1', ...
    'Stretch_4', ...
    'Stretch_3', ...
    'Stretch_2', ...
    'Stretch_1', ...
    'Light_4', ...
    'Light_3', ...
    'Light_2', ...
    'Light_1', ...
    'Kiss_4', ...
    'Kiss_3', ...
    'Kiss_2', ...
    'Kiss_1', ...
    'Kick_4', ...
    'Kick_3', ...
    'Kick_2', ...
    'Kick_1', ...
    'Hit_4', ...
    'Hit_3', ...
    'Hit_2', ...
    'Hit_1', ...
    'Crush_4', ...
    'Crush_3', ...
    'Crush_2', ...
    'Crush_1', ...
    'Console_4', ...
    'Console_3', ...
    'Console_2', ...
    'Console_1'};
userOptions.alternativeConditionLabels=userOptions.conditionLabels;
[userOptions.alternativeConditionLabels{1:64}] = deal(' ');

userOptions.useAlternativeConditionLabels = true;

% What colours should be given to the conditions?
% Verb
userOptions.conditionColours = [repmat([1 0 0], 4,1); ...
    repmat([0 1 0], 4,1); ...
    repmat([0 0 1], 4,1); ...
    repmat([1 1 0], 4,1); ...
    repmat([1 0 1], 4,1); ...
    repmat([0 1 1], 4,1); ...
    repmat([0 0 0], 4,1); ...
    repmat([1 1 1], 4,1) ...
    ];

% Which distance measure to use when calculating first-order RDMs.
userOptions.distance = 'Correlation';

%% Second-order analysis

% Which similarity-measure is used for the second-order comparison.
userOptions.distanceMeasure = 'Spearman';

% How many permutations should be used to test the significance of the fits?  (10,000 highly recommended.)
userOptions.significanceTestPermutations = 1000;
userOptions.nRandomisations=1000;
% Bootstrap options
userOptions.nResamplings = 1000;
userOptions.resampleSubjects = true;
userOptions.resampleConditions = false;
userOptions.plotpValues='=';

% Should RDMs' entries be rank transformed into [0,1] before they're displayed?
userOptions.rankTransform = true;

% Should rubber bands be shown on the MDS plot?
userOptions.rubberbands = true;

% What criterion shoud be minimised in MDS display?
userOptions.criterion = 'metricstress';

% What is the colourscheme for the RDMs?
userOptions.colourScheme = bone(128);

% How should figures be outputted?
userOptions.displayFigures = true;
userOptions.saveFiguresPDF = true;
userOptions.saveFiguresFig = true;
userOptions.saveFiguresPS = false;
% Which dots per inch resolution do we output?
userOptions.dpi = 300;
% Remove whitespace from PDF/PS files?
% Bad if you just want to print the figures since they'll
% no longer be A4 size, good if you want to put the figure
% in a manuscript or presentation.
userOptions.tightInset = true;

% userOptions.forcePromptReply = 's';
