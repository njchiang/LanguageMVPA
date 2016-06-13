function [resp] = slidingBar(window, squaresize)
% this function draws a sliding bar on to PTB Screen window of length
% squaresize

squaresize = 100;

ratenum = [-1 1 ];  %no, yes
rating = round(rand(1,1)*100)/100; % start at random position
ratinginit = rating; % save initial position
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
