% In case something is still ... Remaining from before (which would cause
% our gif not to be correctly displayed... é_è):
clear,clf

% File to open:
% Input Parameters
nPosi = 1;
typeSolver = 'rk4';
rows = 1000; % Number of time steps
nbPlanets = 2; % Number of Planets
printingSteps = 5; % Printing every x time steps ... Don't use it now, but ... Just in case ?
bWantAGif = true; % A Gif or a JPEG ?
% Then processing with the opening
fileToOpen_X = strcat('pos',int2str(nPosi-1),'_',typeSolver,'.dat');
fileToOpen_Y = strcat('pos',int2str(nPosi),'_',typeSolver,'.dat');
fidPosi_X = fopen(fileToOpen_X );
fidPosi_Y = fopen(fileToOpen_Y );
Posi_X = fscanf(fidPosi_X,'%g',[nbPlanets rows]).';
Posi_Y = fscanf(fidPosi_Y,'%g',[nbPlanets rows]).';

if bWantAGif == true 
     % Selection of the name of the file:
    filename = strcat('plot2D_for_' ,typeSolver,'.gif');
    figure
    % Pour tous les X_i
    for i = 1:rows
        for nbPlan = 1:nbPlanets
             if i == 1 ||rem(i,printingSteps) == 0  
                color = rand(1,3);
                if nbPlanets > 1 
                    color = [nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1)];
                end
                plot(Posi_X(i,nbPlan),Posi_Y(i,nbPlan),'color',color,'marker','+');

                hold on % Do we need it ? Don't think so but let's try. Yes we do!
                %drawnow
                frame = getframe(1); % 1 referred to "figure", called before.
                im = frame2im(frame);
                [A,map] = rgb2ind(im,256); % To avoid 3D pictures... For now
               if i == 1 
                   imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',10);
               else
                   imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',10);
               end
                    
            end
        end
    end
else   
    % Do we need a position picture for every planet ? Not really ~
        filename = strcat('plot2D_for_',typeSolver,'.jpeg');
        for nbPlan=1:nbPlanets  % For each planet
            color = rand(1,3);
            if nbPlanets > 1 
                color = [nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1)];               
            end
            plot(Posi_X(:,nbPlan),Posi_Y(:,nbPlan),'color',rand(1,3))
            hold on
        end
        %drawnow % is just a function that allows the drawing of the functions during the process.
        frame = getframe(1);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256); % To avoid 3D pictures
        imwrite(A,map,filename,'jpeg'); 
    
end


fclose(fidPosi_X);
fclose(fidPosi_Y);