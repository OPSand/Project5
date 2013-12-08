fileToOpen = strcat('nbound_0_leapfrog','.dat');
fidPosi = fopen(fileToOpen);
nbSteps = 4;
Plot_File = fscanf(fidPosi,'%g',[10 nbSteps]).';
fclose(fidPosi);
nbParticles = 300;
bWantToSaveJPEG = true;
%filename = strcat('plot_for_',int2str(nbParticles),'_Particles','.jpeg');
filename = 'plotLeapfrog.dat';
% Init of our time vector
t = [1 (nbSteps- 1)];
for i = 1:nbSteps
    t(i)= Plot_File(i,1);
end

% Init of the vector holding the number of bound particles
nbBoundsEk = [0 (nbSteps- 1)];
for i = 1: nbSteps
    nbBoundsEk(i) = Plot_File(i,8);
end
nbBoundsEp = [0 (nbSteps- 1)];
for i = 1: nbSteps
    nbBoundsEp(i) = Plot_File(i,9);
end

nbBoundsVirial = [0 (nbSteps- 1)];
for i = 1: nbSteps
    nbBoundsVirial(i) = nbBoundsEp(i)/nbBoundsEk(i);
end

% And then we plot !
figure(1);
plot(t(:),nbBoundsVirial(:),'color',rand(1,3));
%title(['Plot for ' int2str(nbParticles) ' Particles']);
%xlabel('Time (in t_c_r_u_n_c_h)');
%ylabel('Number of bound particles');
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end
    