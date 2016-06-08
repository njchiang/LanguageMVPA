% PTB Script that runs a full experiment testing the Simon Effect

sentences={'a',


for triali=1:testtrialnum    

    % redraw the backgroup
	Screen('FillRect',window,  bgcolor,winRect);
    indextest1 = testsamples(testidrand(triali));
    rectsizey=sqix(totalsampbins(indextest1,1));%sqxsample(indexrand(triali));
    circlecolori=icirclum(totalsampbins(indextest1,2));
    rectcolor=(icolor(:,totalsampbins(indextest1,3))*255)';%[1 1 1]*128;%((triali-1)*samplenum+(i-1)*30);

    % redraw the frames
    for i=3:3
        squarerange(i,:)=[centeri(i,1)-framesize/2 centeri(i,2)-framesize/2 centeri(i,1)+framesize/2 centeri(i,2)+framesize/2];
        Screen('FrameRect',window,  [255 255 255],squarerange(i,:),[3]);

    end;
            Screen('Flip',window);
    WaitSecs(1);
    
    % draw the random sample
    for i=3:3
        rectpos(i,:) = [centeri(i,1)-rectsizey/2 centeri(i,2)-rectsizey/2 centeri(i,1)+rectsizey/2 centeri(i,2)+rectsizey/2];
	end;
    
    for i=3:3
        % draw the square
        Screen('FillRect', window, rectcolor,rectpos(i,:));        
        % draw circle
        Screen('FillOval',window,  circlecolori,circlepos(i,:));

    end;   
    %%%%%%%%%%%%%%%%%%%%%%%
    
    ratenum = [-1 1 ];  %no, yes
    rating = round(rand(1,1)*100)/100;
    ratinginit = rating;
	ypos = centeri(i,2)+framesize/2+100;
	squaresize = 100;
    xstart = centeri(i,1);%+ratenum(ratei)*(squaresize*1.5);
    squarerange=[xstart-squaresize ypos-squaresize/5 xstart+squaresize ypos+squaresize/5];
    Screen('FrameRect',window,[200 200 200], squarerange,[3]);
    squarerangersp=[xstart-squaresize ypos-squaresize/5 xstart-squaresize+2*squaresize*rating ypos+squaresize/5];
    Screen('FillRect',window,[200 200 200], squarerangersp);
    
    Screen('TextSize',window,40);
    Screen('DrawText',window,'NO(0%)',xstart-squaresize-180, ypos+squaresize/5,[200 200 200]);%
    Screen('DrawText',window,'YES(100%)',xstart+squaresize+30, ypos+squaresize/5,[200 200 200]);%

    ShowCursor(0);
    
    FlushEvents('keyDown');
    
    % Wait for a click and hide the cursor
	Screen('TextSize',window,30);
	Screen('DrawText',window,'Drag mouse (i.e. hold button down) to response',xstart-400,ypos-500,[200 200 200]);
	Screen('DrawText',window,[sprintf('%d',round(rating*100)) '%'],xstart-30, ypos+80,[200 200 200]);%
        Screen('Flip',window);

    while (1)
		[x,y,buttons] = GetMouse;
		if buttons(1)
		  break;
		end
	end
    
    % Loop and track the mouse, drawing the contour
	[theX,theY] = GetMouse;
	thePoints = [theX theY];
    squarerangersp=[xstart-squaresize ypos-squaresize/5 theX ypos+squaresize/5];
    Screen('FillRect',window,[200 200 200], squarerangersp);
    
% 	Screen(theWindow,'DrawLine',255,theX,theY,theX,theY);
	sampleTime = 0.01;
	startTime = GetSecs;
	nextTime = startTime+sampleTime;

	while (1)
		[x,y,buttons] = GetMouse;	
        if ~buttons(1)
            break;
        end
		if (x ~= theX | y ~= theY) & x>=xstart-squaresize & x<=xstart+squaresize
		
            squarerangersp=[xstart-squaresize ypos-squaresize/5 theX ypos+squaresize/5];
            Screen('FillRect',window,[0 0 0], squarerange);
            Screen('DrawText',window,[sprintf('%d',round(rating*100)) '%'],xstart-30, ypos+80,[0 0 0]);%
            Screen('FrameRect',window,[200 200 200], squarerange,[3]);
            Screen('FillRect',window,[200 200 200], squarerangersp);
			theX = x; theY = y;
            
            rating = (theX-(xstart-squaresize))/2/squaresize;
            Screen('DrawText',window,[sprintf('%d',round(rating*100)) '%'],xstart-30, ypos+80,[200 200 200]);%
            
            %%%%%%%%%%%%%%%%%%%
                Screen('Flip',window); % deletes the rest of the screen
            %%%%%%%%%%%%%%%%%%%%
           
		end
		if (GetSecs > nextTime)
			thePoints = [thePoints ; x y];
			nextTime = nextTime+sampleTime;
		end
	end;
    
    rsttest(triali)=rating;
    datarst(triali,1)=rating;
    datarst(triali,2)=totalsampbins(indextest1,1);
    datarst(triali,3)=totalsampbins(indextest1,2);
    datarst(triali,4)=totalsampbins(indextest1,3);
    fprintf(fp,'%d\t%f\t%d\t%d\t%d\t%f\n ',triali,rsttest(triali),totalsampbins(indextest1,1),totalsampbins(indextest1,2),...
    totalsampbins(indextest1,3),ratinginit) ;

    HideCursor;
    SetMouse(center(1),center(2)+200,window);
    WaitSecs(1);    
end;

%{

load('Sentences.mat')
DEBUG=true;
nBlocks=8;
%}


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