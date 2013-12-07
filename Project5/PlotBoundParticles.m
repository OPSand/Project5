fileToOpen = strcat('nbound_0_leapfrog','.dat');
fidPosi = fopen(fileToOpen);
nbSteps = 100;
Plot_File = fscanf(fidPosi,'%g',[2 nbSteps]).';
fclose(fidPosi);
nbParticles = 200;
bWantToSaveJPEG = true;
filename = strcat('plot_for_',int2str(nbParticles),'_Particles','.jpeg');

% Init of our time vector
t = [1 (nbSteps- 1)];
for i = 1:nbSteps
    t(i)= Plot_File(i,1);
end

% Init of the vector holding the number of bound particles
nbBounds = [0 (nbSteps- 1)];
for i = 1: nbSteps
    nbBounds(i) = Plot_File(i,2);
end

% And then we plot !
figure(1);
plot(t(:),nbBounds(:),'color',rand(1,3));
title(['Plot for ' int2str(nbParticles) ' Particles']);
xlabel('Time (in t_c_r_u_n_c_h)');
ylabel('Number of bound particles');
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end
    