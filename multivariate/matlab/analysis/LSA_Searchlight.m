%% RSASearchlight
% simulates fMRI data for a number of subjects and runs searchlight
% analysis using RSA, computes the similarity maps for each subject
% and does group level inference.

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
clear;clc
returnHere = pwd; % We'll come back here later
toolboxRoot = ['Z:/Box Sync/UCLA/Research/LanguageMVPA/code']; addpath(genpath(toolboxRoot));
% Generate a userOptions structure
userOptions = defineUserOptions_LSA(); %edit this
% Generate a simulationOptions structure.
% simulationOptions = simulationOptions_demo_SL();

%% config
searchlightOptions.monitor = false;
searchlightOptions.fisher = true;

searchlightOptions.nConditions=size(betaCorrespondence_LSA(),2);
searchlightOptions.nSessions=size(betaCorrespondence_LSA(),1);
Nsubjects = length(userOptions.subjectNames);

mapsFilename = [userOptions.analysisName, '_fMRISearchlight_Maps.mat'];
RDMsFilename = [userOptions.analysisName, '_fMRISearchlight_RDMs.mat'];
DetailsFilename = [userOptions.analysisName, '_fMRISearchlight_Details.mat'];


%initialize data
userOptions.maskNames={'grayMatter'};
% fullBrainVols = fMRIDataPreparation(betaCorrespondence_Semantics(), userOptions);
% binaryMasks_nS = fMRIMaskPreparation(userOptions);
load('D:\fmri\LanguageMVPA\ImageData\LSA_Masks.mat')
load('D:\fmri\LanguageMVPA\ImageData\cope_ImageData.mat')

% models = constructModelRDMs(modelRDMs_LSA, userOptions);
models = constructModelRDMs(modelRDMs_SL, userOptions);
% models = constructModelRDMs(modelRDMs_ROI_layer, userOptions);
figureRDMs(models, userOptions)
%%compute the correlation maps per subject
% add mask loop?
maskName = 'grayMatter';
userOptions.analysisName='LSA';
for subI = 1:Nsubjects % can parallelize
    subject=userOptions.subjectNames{subI};
    fprintf(['extracting fullBrain volumes for subject %d \n'],subI)
    singleSubjectVols=fullBrainVols.(subject);
    mask=binaryMasks_nS.(subject).(maskName);
    userOptions.searchlightRadius = 12;
    fprintf(['computing correlation maps for subject %d \n'],subI)
    tic
    try
        disp('Parallelizing...')
        [rs, ps, ns, searchlightRDMs.(subject)] = parSearchlightMapping_fMRI(singleSubjectVols, models, mask, userOptions, searchlightOptions);
    catch
        disp('Not parallelizing...')
        [rs, ps, ns, searchlightRDMs.(subject)] = searchlightMapping_fMRI(singleSubjectVols, models, mask, userOptions, searchlightOptions);
    end
    toc
    % NEED TO DELETE STRUCTURALSPATH BEFORE RUNNING BECAUSE I DON'T KNOW HOW IT
    % WORKS
    %     [rs, ps, ns, searchlightRDMs] = fMRISearchlight_jeff(singleSubjectVols, models, mask, userOptions, searchlightOptions);
    for modelI=1:length(models)
        modelName=models(modelI).name;
        gotoDir(userOptions.rootPath, 'Maps/LSA');
       
            fName= strcat('LSA/', subject,  '_', maskName, '_', modelName, '_rMap');
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
    save(['RDMs_', subject, '.mat'], 'searchlightRDMs', '-v7.3');
    clear rs searchlightRDMs
    cd(returnHere)
end %subjectloop
delete(gcp)
%%
%display ROI
% relRoi = sphericalRelativeRoi(userOptions.searchlightRadius,userOptions.voxelSize);
% nVox_searchlight = size(relRoi,1);
% showVoxObj(relRoi+repmat(simulationOptions.effectCen,[nVox_searchlight,1]),1,[3 2 5])
% title(['\bf searchlight with ',num2str(nVox_searchlight),' voxels'])
% handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'SLsimulationSettings'],userOptions);

%% load the previously computed rMaps and concatenate across subjects
if strcmp(userOptions.postprocess, 'MATLAB')
    disp('Postprocessing the subjects. Do you have RAM?')
    % prepare the rMaps:
    for subjectI = 1:Nsubjects
        subData=load([userOptions.rootPath,filesep,'Maps',filesep,'rs_subject',num2str(subjectI),'.mat'])
        rMaps{subjectI} = subData.rs;
        fprintf(['loading the correlation maps for subject %d \n'],subjectI);
    end
    % concatenate across subjects
    for modelI = 1:numel(models)
        for subI = 1:Nsubjects
            thisRs = rMaps{subI};
            thisModelSims(:,:,:,subI) = thisRs(:,:,:,modelI);
        end
        % obtain a pMaps from applying a 1-sided signrank test and also t-test to
        % the model similarities:
        for x=1:size(thisModelSims,1)
            for y=1:size(thisModelSims,2)
                for z=1:size(thisModelSims,3)
                    if mask(x,y,z) == 1
                        [h p1(x,y,z)] = ttest(squeeze(thisModelSims(x,y,z,:)),0,0.05,'right');
                        [p2(x,y,z)] = signrank_onesided(squeeze(thisModelSims(x,y,z,:)));
                    else
                        p1(x,y,z) = NaN;
                        p2(x,y,z) = NaN;
                    end
                end
            end
            disp(x)
        end
        % apply FDR correction
        pThrsh_t = FDRthreshold(p1,0.05,mask);
        pThrsh_sr = FDRthreshold(p2,0.05,mask);
        p_bnf = 0.05/sum(mask(:));
        % mark the suprathreshold voxels in yellow
        supraThreshMarked_t = zeros(size(p1));
        supraThreshMarked_t(p1 <= pThrsh_t) = 1;
        supraThreshMarked_sr = zeros(size(p2));
        supraThreshMarked_sr(p2 <= pThrsh_sr) = 1;  
    end
elseif strcmp(userOptions.postprocess, 'FSL')
    gotoDir(userOptions.rootPath, 'Maps');
    FSL_postprocess(pwd, [userOptions.rootPath '/registration', [userOptions.rootPath '/fnirt'] , [userOptions.rootPath '/masks'], maskName, models, userOptions);
%   FSL_postprocess(Maps_path, warp_path, ref_path, mask_path, maskName, models, userOptions)
    cd(returnHere);
else
    disp('NOT postprocessing')
end


cd(returnHere);