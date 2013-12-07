fileToOpen = strcat('pos',int2str(nPosi),'_',typeSolver,'.dat');
Plot_File = fscanf(fidPosi,'%g',[nbPlanets rows]).';
fclose(fidPosi);
nbSteps = 100;
nbParticles = 100;
bWantToSaveJPEG = false;
filename = strcat('plot_for_',int2str(nbParticles),'_Particles','.jpeg');

% Init of our time vector
t = [0 (nbSteps- 1)];
t(:)= Plot_File(:,1);

% Init of the vector holding the number of bound particles
nbBounds = [0 (nbSteps- 1)];
nbBounds(:) = Plot_File(:,2);

% And then we plot !
figure(1);
plot(t(:),nbBounds(:),'color',rand(1,3));
title(['Plot for ' int2str(nbParticles) ' Particles']);
xlabel('Time');
ylabel('Number of bound particles');
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end
    