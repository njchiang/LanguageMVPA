% PTB Script that runs a full experiment testing the LMVPA

%% parameters
DEBUG=true;
screenRatio = [.1 .1 .8 .85]; %figure this out
%the total number of sentences
nTrialTypes = 4;
%set the indices
trialIdx=[1:nTrialTypes];
%the sentence block
nBlocks=4;
%the random time between showing sentences
displaytime=1:4;
timeout = 2;


%% Part One: Initialization
Screen('Preference', 'VisualDebugLevel', 1);
sinit = input('Subject''s initials: ','s');
try
    seed = str2num(input('Subject seed: ', 's'));
catch
    seed = randi(100);
    fprintf(['Invalid seed, using ' num2str(seed)]);
end

load('sentences.mat');
results = zeros([nTrialTypes, nTrialTypes, 2]);
rng(seed)

%pair the indices and delete the self-repeating pairs
randPairs=CombVec(trialIdx, trialIdx);
randPairs(:,1:nTrialTypes+1:end)=[];

% shuffle the indices here
pairIdx1=Shuffle(1:size(randPairs,1));
pairIdx2=Shuffle(1:size(randPairs,2));
% shuffle the stimuli
shuffledPairs = randPairs(pairIdx1, pairIdx2);

%show the first sentence
%fprintf('S1\n');
input('Press <enter> to begin');
outfilename = ['lmvpa_' sinit];
% [win,Scrnsize] = Screen('OpenWindow',0);
sizeScreen=get(0,'ScreenSize'); 
sizeScreen=[sizeScreen(3:end),sizeScreen(3:end)];
[win, Scrnsize]=Screen('OpenWindow',1,...
    [255 255 255],sizeScreen.*screenRatio);
halfFlip = Screen('GetFlipInterval', win)/2;
KbName('UnifyKeyNames');
centerX = Scrnsize(3)/2;
centerY = Scrnsize(4)/2;
Screen('TextFont',win, 'Arial');
Screen('TextSize',win, 72);

%get the size of the fixation cross
[bounds] = Screen('TextBounds', win, '+');
fixSizeX = bounds(3)/2;

%get the size of the stimuli
[bounds] = Screen('TextBounds', win, 'S2');
stimSizeX = bounds(3)/2;
leftStimX = centerX-200-stimSizeX;
rightStimX = centerX+200-stimSizeX;
if (~DEBUG), ListenChar(2); end

%% Part Two: Data Collection
% for the general condition
if (~DEBUG), HideCursor; end
for blocknumber = 1:nBlocks
    blockLength = size(shuffledPairs,2)/nBlocks;
    theseTrials = shuffledPairs(:, ...
        blockLength*(blocknumber-1) + 1:blockLength); % the trials for this block
    for trialnumber = 1:size(theseTrials, 2)
        displaytime=Shuffle(displaytime);
         % Type of results L/R x on left right 
        WaitSecs(2);   %pause at start of trial, then show fixation

        %%%%%%%%%%% fixSizeX should be sentence 1 offset %%%%%%%%%%%%%%
        Screen('DrawText', win , Sentence{shuffledPairs(1, trialnumber)},...
            centerX-fixSizeX , centerY,[0,0,0]);
        onsetTime1 = Screen('Flip',win);

        %Draw the fixation cross
        Screen('DrawText', win , '+',centerX-fixSizeX ,...
               centerY,[0,0,0]);

        %a flip and show sentence two
        onsetTime2 = Screen('Flip',win,onsetTime1 + displaytime(1)/2 ...
                - halfFlip);
        %%%%%%%%%%% fixSizeX should be sentence 2 offset %%%%%%%%%%%%%%%
        Screen('DrawText', win , Sentence{shuffledPairs(2, trialnumber)}, ...
            centerX-fixSizeX, centerY,[0,0,0]);
        %a flip and show the cross
        onsetTime3 = Screen('Flip',win,onsetTime2 + displaytime(2)/2 ...
                - halfFlip);
        Screen('DrawText', win , '+',centerX-fixSizeX ,...
                centerY,[0,0,0]);
        %a flip and the next round
        onsetTime4 = Screen('Flip',win,onsetTime3 + displaytime(3)/2 ...
                - halfFlip);
        Screen('DrawText', win , 'Question goes here', ...
            centerX-fixSizeX, centerY,[0,0,0]);
        Starttime  = Screen('Flip',win,onsetTime4 + displaytime(4)/2 ...
                - halfFlip);
        Nowtime = Starttime;
        % collect rating here:
        % [rating, respTime] = slidingBar(window); % this won't work
        rating = rand;
        respTime = rand;
        results(shuffledPairs(1, trialnumber), shuffledPairs(2, trialnumber), 1) = rating;
        results(shuffledPairs(1, trialnumber), shuffledPairs(2, trialnumber), 2) = respTime;

    end
 Screen('DrawText', win , 'End of block... press Enter to continue', ...
            centerX-fixSizeX, centerY,[0,0,0]);
 Starttime  = Screen('Flip',win,onsetTime4 + displaytime(4)/2 ...
        - halfFlip); 
% get some sort of button press here
input('Press <enter>');

end
%% Part Three: Cleanup and File save
ShowCursor
if (~DEBUG), ListenChar(0); end
sca;
save(outfilename, 'results', 'shuffledPairs', 'seed', 'sinit');