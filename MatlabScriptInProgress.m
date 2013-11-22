% File to open:
nPosi = 0;
typeSolver = 'rk4';
fileToOpen = strcat('pos',int2str(nPosi),'_',typeSolver,'.dat');
fidPosi = fopen(fileToOpen);
rows = 1000; % Number of time steps
nbPlanets = 2; % Number of Planets
printingSteps = 50; % Printing every x time steps ... Don't use it now, but ... Just in case ?
bWantAGif = true; % A Gif or a JPEG ?
Posi = fscanf(fidPosi,'%g',[nbPlanets rows]).';

% Init of our time vector
t = [0 (rows - 1)];
for i= 1:rows
    t(i)= i*1/rows;
end

machin = 'Truc';

if bWantAGif == true  
    % Selection of the name of the file:
    filename = strcat('plot_for_' ,typeSolver,'.gif');
    figure
    % For every time step
        for time = 1:rows   
            for nbPlan = 1:nbPlanets
                if time == 1
                      imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',30);
                else
                    if rem(time,printingSteps) == 0
                        color = rand(1,3);
                        if nbPlanets > 1 
                            color = [nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1)];
                        end
                        plot(t(time),Posi(time,nbPlan),'color',color,'marker','+');

                        hold on % Do we need it ? Don't think so but let's try. Yes we do!
                        drawnow
                        frame = getframe(1); % 1 referred to "figure", called before.
                        im = frame2im(frame);
                        [A,map] = rgb2ind(im,256); % To avoid 3D pictures... For now
                       
                    end
                end
            end
            % Then let's turn everything into a gif
%             if time == 1 % After the first step, we create the gif
%                 imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',30);
%             else % Then, if we are in our printing step
                if rem(time,printingSteps) == 0 % We append the new gif to the current one
                     imwrite(A,map,filename,'gif','WriteMode','append')%,'DelayTime',10);
                 end
%             end     
        end
else   
    % Do we need a position picture for every planet ? Not really ~
        filename = strcat('plot_for_',typeSolver,'.jpeg');
        for nbPlan=1:nbPlanets  % For each planet
            color = rand(1,3);
            if nbPlanets > 1 
                color = [nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1)];               
            end
            plot(t,Posi(:,nbPlan),'color',color) % We plot the current position
            hold on
        end
        %drawnow % is just a function that allows the drawing of the functions during the process.
        frame = getframe(1);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256); % To avoid 3D pictures
        imwrite(A,map,filename,'jpeg'); 
    
end

fclose(fidPosi);


%     for nbPlan = 1:nbPlanets
%     filename = strcat('plot',int2str(nbPlan),'for' ,typeSolver,'.gif');
%     %figure
%     % For every time step
%         for j = 1:rows                  
%             plot(t,Posi(j,nbPlan),'color',rand(1,3))
%             drawnow
%         end
%         frame = getframe(1);
%             im = frame2im(frame);
%             [A,map] = rgb2ind(im,256); % To avoid 3D pictures
%             if nbPlan == 1 % After the first step, we create the gif
%                 imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%             else % Then, if we are in our printing step
%                 %if rem(j,printingSteps) == 0 % We append the new gif to the current one
%                 imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
%                 %end
%             end
%         
%         
%     end

% Now I want to hold on every ... Plot_xxxY.dat and draw everything on the
% same file.
% We have a posx File for every Dimension. So let's focus on the plot1 for
% now.
% name = '';
% % Plotting and printing
% for i=1:nbPlanets
%  % plot = plot(X(:,i),Y(:,i),'color',rand(1,3));
%   
%   for j = 2:rows  
%        %plot(X_i(j),t(j));
%     if rem(j,printingSteps) == 0
%             name = strcat('plot',int2str(i),'.jpeg');
%             %print(plot,name);
%     end
%      name; % Just to check if the name is as expected
%   end
% end