% In case something is still ... Remaining from before (which would cause
% our gif not to be correctly displayed... é_è):
clear,clf

% File to open:
% Input Parameters
nPosi = 2; % From 0 to nPosi
typeSolver = 'rk4';
rows = 1000; % Number of time steps
nbPlanets = 100; % Number of Planets
printingSteps = 20; % Printing every x time steps ... Don't use it now, but ... Just in case ?
bWantAGif = true; % A Gif or a JPEG ?
% Then processing with the opening
fileToOpen_X = strcat('pos',int2str(nPosi-2),'_',typeSolver,'.dat');
fileToOpen_Y = strcat('pos',int2str(nPosi-1),'_',typeSolver,'.dat');
fileToOpen_Z = strcat('pos',int2str(nPosi),'_',typeSolver,'.dat');
fidPosi_X = fopen(fileToOpen_X );
fidPosi_Y = fopen(fileToOpen_Y );
fidPosi_Z = fopen(fileToOpen_Z );
Posi_X = fscanf(fidPosi_X,'%g',[nbPlanets rows]).'
Posi_Y = fscanf(fidPosi_Y,'%g',[nbPlanets rows]).'
Posi_Z = fscanf(fidPosi_Z,'%g',[nbPlanets rows]).'

% And finally closing the handle on the open documents
fclose(fidPosi_X);
fclose(fidPosi_Y);
fclose(fidPosi_Z);

if bWantAGif == true
      % Selection of the name of the file:
    filename = strcat('plot3D_for_' ,typeSolver,'.gif');
    figure(1)
    grid on;
h = plot3(Posi_X(:,1),Posi_Y(:,1),Posi_Z(:,1),'r.');

xlim([-40 40]);
ylim([-40 40]);
zlim([-40 40]);
for i = 2: rows
    set(h,'XData',Posi_X(i,:),'YData',Posi_Y(i,:),'ZData',Posi_Z(i,:));
    
    grid on;
    drawnow
    frame = getframe(1);
                im = frame2im(frame);
                [A,map] = rgb2ind(im,256); % To avoid 3D pictures... For now
                
         if i == 2
           imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
         else
            imwrite(A,map,filename,'gif','WriteMode','append');
        end
    pause(0.1);
end
end


%     % Pour tous les X_i
%     for i = 2:rows
%              xlim([-100 100]);
%                 ylim([-100 100]);
%                 zlim([-100 100]);
%             color = rand(1,3);
%             if i == 1 ||rem(i,printingSteps) == 0  
%                 for nbPlan = 1:nbPlanets    
%                 if nbPlanets > 1 
%                     color = [nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1)];
%                 end
%                 
%                 plot3(Posi_X(i,nbPlan),Posi_Y(i,nbPlan),Posi_Z(i,nbPlan),'color',color,'marker','o');

%                 hold on % Do we need it ? Don't think so but let's try. Yes we do!
%                 drawnow
%                 end
%                 title(['Plot for ', typeSolver]);
%                 xlabel('X');
%                 ylabel('Y');
%                 zlabel('Z');
%                 grid on;
%                 %frame = getframe(1);
%                 %im = frame2im(frame);
%                 %[A,map] = rgb2ind(im,256); % To avoid 3D pictures... For now
%                 %if i == 1 
%                 %    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
%                 %else
%                  %   imwrite(A,map,filename,'gif','LoopCount',Inf,'WriteMode','overwrite');
%                     %imwrite(A,map,filename,'gif','WriteMode','append');
%                 %end
%                     
%             
%             %frame(i) = getframe(1); % 1 referred to "figure", called before.s
%             end
%         hold off;
%     end
%     %movie(frame,10);
%     close all
% else
%     % Do we need a position picture for every planet ? Not really ~
%     filename = strcat('plot3D_for_',typeSolver,'.jpeg');
%     for nbPlan=1:nbPlanets  % For each planet
%         color = rand(1,3);
%         if nbPlanets > 1 
%             color = [nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1) nbPlan/(nbPlanets+1)];               
%         end
%         plot3(Posi_X(:,nbPlan),Posi_Y(:,nbPlan),Posi_Z(:,nbPlan),'color',rand(1,3))
%         hold on
%     end
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     title(['Plot for ', typeSolver]);
%     %drawnow % is just a function that allows the drawing of the functions during the process.
%     frame = getframe(1);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256); % To avoid 3D pictures
%     imwrite(A,map,filename,'jpeg'); 
% end