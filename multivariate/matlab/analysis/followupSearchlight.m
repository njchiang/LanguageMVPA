%% RSASearchlight
% simulates fMRI data for a number of subjects and runs searchlight
% analysis using RSA, computes the similarity maps for each subject
% and does group level inference.

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
clear;clc
% cd ..
toolboxRoot = ['~/Box Sync/UCLA/Research/LanguageMVPA/code']; addpath(genpath(toolboxRoot));
% Generate a userOptions structure
cd /Volumes/fmri/LanguageMVPA
returnHere = pwd; % We'll come back here later

userOptions = defineUserOptions_LSA(); %edit this
userOptions.searchlightRadius=6;
userOptions.analysisName='Across_Searchlight';
% userOptions.rootPath = [pwd,filesep];
% userOptions.analysisName = 'Searchlight';

% Generate a simulationOptions structure.
% simulationOptions = simulationOptions_demo_SL();

%% config
searchlightOptions.monitor = false;
searchlightOptions.fisher = true;

searchlightOptions.nConditions=size(betaCorrespondence_Semantics(),2);
searchlightOptions.nSessions=size(betaCorrespondence_Semantics(),1);
Nsubjects = length(userOptions.subjectNames);

mapsFilename = [userOptions.analysisName, '_fMRISearchlight_Maps.mat'];
RDMsFilename = [userOptions.analysisName, '_fMRISearchlight_RDMs.mat'];
DetailsFilename = [userOptions.analysisName, '_fMRISearchlight_Details.mat'];


%% loading structurals: I DON'T DO THIS
% load([returnHere,filesep,'sampleMask_org.mat'])
% load([returnHere,filesep,'anatomy.mat']);% load the resliced structural image
% warpFlags.interp = 1;
% warpFlags.wrap = [0 0 0];
% userOptions.voxelSize = [3 3 3];
% warpFlags.vox = userOptions.voxelSize; % [3 3 3.75]
% warpFlags.bb = [-78 -112 -50; 78 76 85];
% warpFlags.preserve = 0;

%initialize data
userOptions.maskNames={'grayMatter'};
% binaryMasks_nS = fMRIMaskPreparation(userOptions);
load('ImageData/LSA_Masks.mat')
load('ImageData/LSA_ImageData.mat')
% load('E:\LanguageMVPA\ImageData\Sub_Semantics_ImageData.mat')
% userOptions.maskPath='E:\LanguageMVPA\masks\[[subName]]_[[maskName]].img';
% % userOptions.maskNames={'LH_InferiorFrontalGyrus'};
% binaryMasks_nS = fMRIMaskPreparation(userOptions);
% load('ImageData/SearchlightMusic_Masks')
% models = constructModelRDMs(modelRDMs_searchlight2, userOptions);
models(1).name='Lin_L2PAnim';
models(1).label=[1*ones(1,4) 2*ones(1,4) 1*ones(1,4) 2*ones(1,8) 1*ones(1,8) 2*ones(1,4)];
models(2).name='Lin_P2LAnim';
models(2).label=models(1).label;
models(3).name='Lin_L2LAnim';
models(3).label=models(1).label;
models(4).name='Lin_P2PAnim';
models(4).label=models(1).label;
models(5).name='Lin_L2PFam';
models(5).label=[1*ones(1,4) 2*ones(1,8) 1*ones(1,4) 3*ones(1,8) 4*ones(1,8)];
models(6).name='Lin_P2LFam';
models(6).label=models(5).label;
models(7).name='Lin_L2LFam';
models(7).label=models(5).label;
models(8).name='Lin_P2PFam';
models(8).label=models(5).label;
models(9).name='Lin_L2PVerb';
models(9).label=[1*ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4) 6*ones(1,4) 7*ones(1,4) 8*ones(1,4)];
models(10).name='Lin_P2LVerb';
models(10).label=models(9).label;
models(11).name='Lin_L2LVerb';
models(11).label=models(9).label;
models(12).name='Lin_P2PVerb';
models(12).label=models(9).label;
searchlightModels=[models(1).label; models(5).label; models(9).label];

models(13).name='RBF_L2PAnim';
models(14).name='RBF_P2LAnim';
models(15).name='RBF_L2LAnim';
models(16).name='RBF_P2PAnim';
models(17).name='RBF_L2PFam';
models(18).name='RBF_P2LFam';
models(19).name='RBF_L2LFam';
models(20).name='RBF_P2PFam';
models(21).name='RBF_L2PVerb';
models(22).name='RBF_P2LVerb';
models(23).name='RBF_L2LVerb';
models(24).name='RBF_P2PVerb';
%%compute the correlation maps per subject
% add mask loop?
% maskName = userOptions.maskNames;
maskName='grayMatter'; % set the mask
% maskName='LH_InferiorFrontalGyrus'; % set the mask

% parpool open

for subI = 1:Nsubjects % can parallelize
% for subI = 14:Nsubjects % can parallelize
    subject=userOptions.subjectNames{subI};
    fprintf(['extracting fullBrain volumes for subject %d \n'],subI)
    singleSubjectVols=fullBrainVols.(subject);
    
    userOptions.searchlightRadius = 6;
    mask = binaryMasks_nS.(subject).(maskName);
    fprintf(['computing correlation maps for subject %d \n'],subI)
    %     [rs, ps, ns, searchlightRDMs.(subject)] = searchlightMapping_fMRI(singleSubjectVols, models, mask, userOptions, searchlightOptions);
    tic
    [rs, ps, ns] = searchlight_LMVPA_across_SVM(singleSubjectVols, searchlightModels, mask, userOptions, searchlightOptions);
    toc
    delete(gcp)
    % NEED TO DELETE STRUCTURALSPATH BEFORE RUNNING BECAUSE I DON'T KNOW HOW IT
    % WORKS
    %     [rs, ps, ns, searchlightRDMs] = fMRISearchlight_jeff(singleSubjectVols, models, mask, userOptions, searchlightOptions);
    
    for modelI=1:length(models)
        modelName=models(modelI).name;
        gotoDir(userOptions.rootPath, 'Maps/Across');
        
        fName= strcat('Across/', subject,  '_', maskName, '_', modelName, '_rMap');
        if class(fName)=='cell'
            writeOpts.name=fName{1};
        else
            writeOpts.name=fName;
        end
        writeOpts.description=[subject '_R-Map'];
        writeOpts.template=[userOptions.rootPath '/template_brain.hdr'];
        writeMe=rs(:,:,:,modelI);
        write_brainMap(writeMe, userOptions, writeOpts);
    end
    save(['rs_',subject,'.mat'],'rs');
%     save(['RDMs_', subject, '.mat'], 'searchlightRDMs', '-v7.3');
    clear rs 
    cd(returnHere)
    
end %subjectloop

