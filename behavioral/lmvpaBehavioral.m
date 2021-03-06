function [] = lmvpaBehavioral()
% PTB Script that runs a full experiment testing the LMVPA
% to do... only probe upper or lower triangle
% maybe we want to incrementally save after each block...
%% parameters to adjust
DEBUG_MODE=false; % true for half screen)
% number of blocks... probably need to set this per trial.
nBlocks=2;

textSize=32;
screenRatio = [0 0 .5 .5]; % offsetX offsetY sizeX sizeY



% set the indices of sentences to be used
displayTimeMean=2; % randomly sample from exponential with expectation 2s
fixDurMean = 1; % randomly sample from exponential with expectation 1s
barSize = [300 30];
trialSelector = 'upper'; % 'upper', 'lower', or 'all'
sentenceDisplayRange=[1.5 3];
fixDisplayRange=[.5 1.5];

%% Part One: Initialization
% if strcmpi(DEBUG_MODE, 'true')
% Screen('Preference', 'ConserveVRAM', 64);
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference', 'VisualDebugLevel', 1);
% end
sinit = input('Subject''s initials: ','s');

cond = input('Stimulus list: ', 's');
% probably make other matrices for syntax and semantics
load([cond '.mat']);
% the total number of sentences
nTrialTypes = length(Sentence);
outfilename = ['lmvpa_' cond '_' sinit];

try
    seed = str2double(input('Subject seed: ', 's'));
catch
    seed = randi(100);
    fprintf(['Invalid seed, using ' num2str(seed)]);
end

results = zeros([nTrialTypes, nTrialTypes, 2]);
rng(seed)
trialIdx=[1:nTrialTypes];
shuffledPairs = shuffleIdx(trialIdx, nTrialTypes);
% top half or bottom half of matrix depending on sign
if strcmp('upper', trialSelector)
    triangleSelector = shuffledPairs(1,:) > shuffledPairs(2,:);
elseif strcmp('lower', trialSelector)
    triangleSelector = shuffledPairs(1,:) < shuffledPairs(2,:);
else
    triangleSelector = ones(1, size(shuffledPairs,2));
end
shuffledPairs = shuffledPairs(:, triangleSelector);

fixDurVec =  exprnd(fixDurMean, size(shuffledPairs));
fixDurVec(fixDurVec < fixDisplayRange(1)) = fixDisplayRange(1);
fixDurVec(fixDurVec > fixDisplayRange(2)) = fixDisplayRange(2);
displayTimeVec = exprnd(displayTimeMean, size(shuffledPairs));
displayTimeVec(displayTimeVec < sentenceDisplayRange(1)) = sentenceDisplayRange(1);
displayTimeVec(displayTimeVec > sentenceDisplayRange(2)) = sentenceDisplayRange(2);


input('Press <enter> to begin');
fullScreenSize=get(0,'ScreenSize');
HideCursor(0);
if DEBUG_MODE
    
    [win, ScreenSize]=Screen('OpenWindow',0,...
        [255 255 255],fullScreenSize.*screenRatio);
else
    [win, ScreenSize]=Screen('OpenWindow',0,...
        [255 255 255]);
end
KbName('UnifyKeyNames');
Screen('TextFont',win, 'Arial');
Screen('TextSize',win, textSize);

if (~DEBUG_MODE), ListenChar(2); end
% INSTRUCTIONS GO HERE

%% Part Two: Data Collection
% for the general condition
if (~DEBUG_MODE), HideCursor; end
for blocknumber = 1:nBlocks
    blockLength = size(shuffledPairs,2)/nBlocks;
    theseTrials = shuffledPairs(:, ...
        blockLength*(blocknumber-1) + [1:blockLength]); % the trials for this block
    theseDisplayTimes = displayTimeVec(:, ...
        blockLength*(blocknumber-1) + [1:blockLength]);
    theseFixTimes = fixDurVec(:, ...
        blockLength*(blocknumber-1) + [1:blockLength]);
    for trialnumber = 1:size(theseTrials, 2)
        % displaytime=Shuffle(displaytime);
        HideCursor(0);
        trialStartTime = Screen('Flip', win);
        
        %%%%%%%%%%% each trial %%%%%%%%%%%%%%
        % display first sentence
        trial1Time = displayText(win, ScreenSize, ...
            Sentence{shuffledPairs(1, trialnumber)}, ...
            trialStartTime, theseFixTimes(1, trialnumber));
        % first fixation
        fix1Time = displayText(win, ScreenSize, '+', ...
            trial1Time, theseDisplayTimes(1, trialnumber));
        % second sentence
        trial2Time = displayText(win, ScreenSize, ...
            Sentence{shuffledPairs(2, trialnumber)}, ...
            fix1Time, theseFixTimes(2, trialnumber));
        % second fixation
        fix2Time = displayText(win, ScreenSize, '+', ...
            trial2Time, theseDisplayTimes(2, trialnumber));
        % collect rating here:
        nowTime=GetSecs;
        [rating, ~, ratinginit] = slidingBar(win, barSize, ScreenSize);
        respTime = GetSecs-nowTime;
        % reset text parameters
        Screen('TextFont',win, 'Arial');
        Screen('TextSize',win, textSize);
        % save results
        results(shuffledPairs(1, trialnumber), shuffledPairs(2, trialnumber), 1) = rating;
        results(shuffledPairs(1, trialnumber), shuffledPairs(2, trialnumber), 2) = respTime;
        
    end
    endBlockText='End of block... press Enter to continue';
    displayText(win, ScreenSize, endBlockText);
    % incremental save
    save(outfilename, 'results', 'shuffledPairs', 'seed', 'sinit');
    % press a button to advance
    % probably should use Kb...
    %     input('Press <enter>');
    cResponse = 0;
    while (cResponse == 0)
        [keyDown, ~, ~] = KbCheck(-1);
        if(keyDown)
            cResponse = 1;
        end
    end
end
%% Part Three: Cleanup and File save
ShowCursor
if (~DEBUG_MODE), ListenChar(0); end
save(outfilename, 'results', 'shuffledPairs', 'seed', 'sinit');
endExperimText='Thank you for participating. Press <enter> to exit';
displayText(win, ScreenSize, endExperimText);
% use kb here
% input('Press <enter>');

cResponse = 0;
while (cResponse == 0)
    [keyDown, ~, ~] = KbCheck(-1);
    if(keyDown)
        cResponse = 1;
    end
end
sca;

end

function [flipTime] = displayText(window, ScreenSize, text, onsetTime, duration)
halfFlip = Screen('GetFlipInterval', window)/2;
[bounds] = Screen('TextBounds', window, text);
stimSizeX = bounds(3)/2;
stimSizeY = bounds(4)/2;
centerX = ScreenSize(3)/2;
centerY = ScreenSize(4)/2;
Screen('DrawText', window , text, ...
    centerX-stimSizeX, centerY-stimSizeY,[0,0,0]);
if nargin > 3
    flipTime = Screen('Flip',window,onsetTime + duration ...
        - halfFlip);
else
    flipTime = Screen('Flip', window);
end
end

function shuffledPairs = shuffleIdx(trialIdx, nTrialTypes)
%pair the indices and delete the self-repeating pairs
randPairs=CombVec(trialIdx, trialIdx);
randPairs(:,1:nTrialTypes+1:end)=[];
% shuffle the indices here
pairIdx1=Shuffle(1:size(randPairs,1));
pairIdx2=Shuffle(1:size(randPairs,2));
% shuffle the stimuli
shuffledPairs = randPairs(pairIdx1, pairIdx2);
end