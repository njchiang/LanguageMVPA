function [] = lmvpaBehavioral()
% PTB Script that runs a full experiment testing the LMVPA

%% parameters to adjust
DEBUG=true;
textSize=48;
screenRatio = [0 0 .5 .5]; % offsetX offsetY sizeX sizeY
% the total number of sentences
nTrialTypes = 4;
% set the indices of sentences to be used
trialIdx=[1:nTrialTypes];
nBlocks=4;
displaytime=2;
fixdur = 1;
barSize = [100 20];

%% Part One: Initialization
Screen('Preference', 'VisualDebugLevel', 1);
sinit = input('Subject''s initials: ','s');
outfilename = ['lmvpa_' sinit];
try
    seed = str2double(input('Subject seed: ', 's'));
catch
    seed = randi(100);
    fprintf(['Invalid seed, using ' num2str(seed)]);
end

load('sentences.mat');
results = zeros([nTrialTypes, nTrialTypes, 2]);
rng(seed)
shuffledPairs = shuffleIdx(trialIdx, nTrialTypes);

input('Press <enter> to begin');
fullScreenSize=get(0,'ScreenSize');
[win, ScreenSize]=Screen('OpenWindow',0,...
    [255 255 255],fullScreenSize.*screenRatio);
KbName('UnifyKeyNames');
Screen('TextFont',win, 'Arial');
Screen('TextSize',win, textSize);

if (~DEBUG), ListenChar(2); end

%% Part Two: Data Collection
% for the general condition
if (~DEBUG), HideCursor; end
for blocknumber = 1:nBlocks
    blockLength = size(shuffledPairs,2)/nBlocks;
    theseTrials = shuffledPairs(:, ...
        blockLength*(blocknumber-1) + [1:blockLength]); % the trials for this block
    for trialnumber = 1:size(theseTrials, 2)
        % displaytime=Shuffle(displaytime);
        trialStartTime = Screen('Flip', win);
        
        %%%%%%%%%%% each trial %%%%%%%%%%%%%%
        trial1Time = displayText(win, ScreenSize, ...
            Sentence{shuffledPairs(1, trialnumber)}, trialStartTime, fixdur);
        
        fix1Time = displayText(win, ScreenSize, '+', ...
            trial1Time, displaytime);
        
        trial2Time = displayText(win, ScreenSize, ...
            Sentence{shuffledPairs(2, trialnumber)}, fix1Time, fixdur);
        
        fix2Time = displayText(win, ScreenSize, '+', trial2Time, displaytime);
        
        % collect rating here:
        nowTime=GetSecs;
        [rating, ~, ratinginit] = slidingBar(win, barSize, ScreenSize);
        
        respTime = GetSecs-nowTime;
        % reset text parameters
        Screen('TextFont',win, 'Arial');
        Screen('TextSize',win, textSize);
        
        results(shuffledPairs(1, trialnumber), shuffledPairs(2, trialnumber), 1) = rating;
        results(shuffledPairs(1, trialnumber), shuffledPairs(2, trialnumber), 2) = respTime;
        
    end
    endBlockText='End of block... press Enter to continue';
    displayText(win, ScreenSize, endBlockText);
    % press a button to advance
    input('Press <enter>');
    
end
%% Part Three: Cleanup and File save
ShowCursor
if (~DEBUG), ListenChar(0); end
save(outfilename, 'results', 'shuffledPairs', 'seed', 'sinit');
input('Thank you for participating. Press <enter> to exit');
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