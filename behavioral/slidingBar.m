function [rating, nextTime, ratinginit] = slidingBar(window, squaresize, ScreenSize)
% this function draws a sliding bar on to PTB Screen window of length
% squaresize
textSize = 40;
rating = round(rand(1,1)*100)/100; % start at random position
ratinginit = rating; % save initial position
ypos = ScreenSize(4)/2; % random starting point
xstart = ScreenSize(3)/2; % random starting point
squarerange=[xstart-squaresize(1) ypos-squaresize(2) xstart+squaresize(1) ypos+squaresize(2)];
Screen('FrameRect',window,[200 200 200], squarerange,[3]);
squarerangersp=[xstart-squaresize(1) ypos-squaresize(2) xstart-squaresize(1)+2*squaresize(1)*rating ypos+squaresize(2)];
Screen('FillRect',window,[200 200 200], squarerangersp);

Screen('TextSize',window, textSize);
% Screen('DrawText',window,'NO(0%)',xstart-squaresize(1)-180, ypos+squaresize(2),[200 200 200]);%
% Screen('DrawText',window,'YES(100%)',xstart+squaresize(1)+30, ypos+squaresize(2),[200 200 200]);%

ShowCursor(0);

FlushEvents('keyDown');

% Wait for a click and hide the cursor
% Screen('DrawText',window,'Drag mouse (i.e. hold button down) to respond',xstart-400,ypos-500,[200 200 200]);
% Screen('DrawText',window,[sprintf('%d',round(rating*100)) '%'],xstart-30, ypos+80,[200 200 200]);%
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
squarerangersp=[xstart-squaresize(1) ypos-squaresize(2) theX ypos+squaresize(2)];
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
    if (x ~= theX || y ~= theY) && x>=xstart-squaresize(1) && x<=xstart+squaresize(1)
        
        squarerangersp=[xstart-squaresize(1) ypos-squaresize(2) theX ypos+squaresize(2)];
        Screen('FillRect',window,[0 0 0], squarerange);
%         Screen('DrawText',window,[sprintf('%d',round(rating*100)) '%'],xstart-30, ypos+80,[0 0 0]);%
        Screen('FrameRect',window,[200 200 200], squarerange,[3]);
        Screen('FillRect',window,[200 200 200], squarerangersp);
        theX = x; theY = y;
        
        rating = (theX-(xstart-squaresize(1)))/2/squaresize(1);
%         Screen('DrawText',window,[sprintf('%d',round(rating*100)) '%'],xstart-30, ypos+80,[200 200 200]);%
        
        %%%%%%%%%%%%%%%%%%%
        Screen('Flip',window); % deletes the rest of the screen
        %%%%%%%%%%%%%%%%%%%%
        
    end
    if (GetSecs > nextTime)
        thePoints = [thePoints ; x y];
        nextTime = nextTime+sampleTime;
    end
end;

HideCursor(0);

