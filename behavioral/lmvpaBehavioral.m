% PTB Script that runs a full experiment testing the Simon Effect

DEBUG=true;
nBlocks=8;

%% Part One: Initialization
Screen('Preference', 'VisualDebugLevel', 1);
sinit = input('Subject''s initials: ','s');
try
    seed = str2num(input('Subject seed: ', 's'));
catch
    seed = randi(100);
    fprintf(['Invalid seed, using ' num2str(seed)]);
end
% change this to task instructions
fprintf('Fixate on "+", press <A> or <;> for L or R, resp.\n');
input('Press <enter> to begin');
outfilename = ['lmvpa_' sinit];
rng(seed)
% [win,Scrnsize] = Screen('OpenWindow',0);
sizeScreen=get(0,'ScreenSize'); 
sizeScreen=[sizeScreen(3:end),sizeScreen(3:end)];
[win, Scrnsize]=Screen('OpenWindow',0,...
    [255 255 255], sizeScreen.*[.2 .5 .8 .85]);
halfFlip = Screen('GetFlipInterval', win)/2;
KbName('UnifyKeyNames');
centerX = Scrnsize(3)/2;
centerY = Scrnsize(4)/2;
Screen('TextFont',win, 'Arial');
Screen('TextSize',win, 72);
timeout = 2;
% change this to sentences
[results(1:4).side] = deal('L','R','L','R');
[results(1:4).stim] = deal('L','L','R','R');
[results(1:4).comp] = deal('C','I','I','C');
% Prepare data fields for each type of trial.
[results(1:4).RT] = deal([]);
[results(1:4).error] = deal(0);
%get the size of the fixation cross
[bounds] = Screen('TextBounds', win, '+');
fixSizeX = bounds(3)/2;
%get the size of the stimuli
[bounds] = Screen('TextBounds', win, 'L');
stimSizeX = bounds(3)/2;
leftStimX = centerX-200-stimSizeX;
rightStimX = centerX+200-stimSizeX;
if (~DEBUG), ListenChar(2); end

%% Part Two: Data Collection
if (~DEBUG), HideCursor; end
for blocknumber = 1:nBlocks
    for typenum = randperm(4); % Type of results L/R x on left right
        WaitSecs(2);   %pause at start of trial, then show fixation
        Screen('DrawText', win , '+',centerX-fixSizeX ,...
            centerY,[0,0,0]);
        onsetTime1 = Screen('Flip',win);
        %Draw the fixation cross
        Screen('DrawText', win , '+',centerX-fixSizeX ,...
            centerY,[0,0,0]);
        %Show the stimulus on the left or right side
        if results(typenum).side == 'L'
            Screen('DrawText', win , results(typenum).stim,...
                leftStimX,centerY,[0,0,0]);
        else
            Screen('DrawText', win , results(typenum).stim,...
                rightStimX,centerY,[0,0,0]);
        end
        Starttime  = Screen('Flip',win,onsetTime1 + 1.0 ...
            - halfFlip);
        Nowtime = Starttime;
        responseGiven = 0;
        response = 0;
        %collect a response with a timeout
        while ((Nowtime < Starttime + timeout) && (responseGiven == 0))
            %Check for a response
            [keyDown,secs,keyCode]=KbCheck(-1);
            if(keyDown)
                responseGiven = 1;
                response = KbName(keyCode);
            end
            Nowtime = GetSecs;       % check the current time
        end
        thisRT = secs-Starttime;   %compute the reaction time
        
        if (response(1)=='a')   %convert the response into L or R
            thisResp = 'L';
        elseif (response(1) == ';')
            thisResp = 'R';
        else
            thisResp = 'X';
        end
        if (results(typenum).stim == thisResp)
            results(typenum).RT = [results(typenum).RT thisRT];
        else
            results(typenum).error = results(typenum).error + 1;
            Beeper;
        end
        Screen('Flip',win);
    end
end

%% Part Three: Cleanup and File save
ShowCursor
if (~DEBUG), ListenChar(0); end
sca;
save(outfilename,'results');