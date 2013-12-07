% In case something is still ... Remaining from before (which would cause
% our gif not to be correctly displayed... é_è):
clear,clf

% File to open:
% Input Parameters
nPosi = 2; % From 0 to nPosi
typeSolver = 'rk4';
rows =  300*365; % Number of time steps
nbPlanets = 2; % Number of Planets
printingSteps = 10; % Printing every x time steps ... Don't use it now, but ... Just in case ?
bWantAGif = false; % A Gif or a JPEG ?
bWantToSave = false; % This considerably slows down everything, but well ... !
% Then processing with the opening
fileToOpen_X = strcat('pos',int2str(nPosi-2),'_',typeSolver,'.dat');
fileToOpen_Y = strcat('pos',int2str(nPosi-1),'_',typeSolver,'.dat');
fileToOpen_Z = strcat('pos',int2str(nPosi),'_',typeSolver,'.dat');
fidPosi_X = fopen(fileToOpen_X );
fidPosi_Y = fopen(fileToOpen_Y );
fidPosi_Z = fopen(fileToOpen_Z );
Posi_X = fscanf(fidPosi_X,'%g',[nbPlanets rows]).';
Posi_Y = fscanf(fidPosi_Y,'%g',[nbPlanets rows]).';
Posi_Z = fscanf(fidPosi_Z,'%g',[nbPlanets rows]).';

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
    xlim([-2.0e+11 2.0e+11]);
    ylim([-2.0e-04    2.0e-04]);
    zlim([-2.0e+12 2.0e+12]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Plot with ', typeSolver]);
    for i = 2: rows
        %if bWantToSave == true
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
%         else
%             set(h,'XData',Posi_X(i,:),'YData',Posi_Y(i,:),'ZData',Posi_Z(i,:));
%             grid on;
%             drawnow
%         end
        pause(0.1);
    end
else
    figure(1)
    filename = strcat('plot3D_for_',typeSolver,'.jpeg');
    for i = 1: rows
       
            % Do we need a position picture for every planet ? Not really ~
            plot3(Posi_X(i,:),Posi_Y(i,:),Posi_Z(i,:),'r.');
            grid on;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            hold on;
    end
    title(['Plot with ', typeSolver]);
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end