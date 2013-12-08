fileLeap = strcat('plotLeapfrog.dat');
fidLeap = fopen(fileLeap);
Plot_File = fscanf(fidLeap,'%g',[18 4]).';
fclose(fidLeap);

simID = 3;
nbParticles = 2000;
nbSteps = 67;
fileRad = strcat('radial_after_',num2str(simID),'_leapfrog.dat');
fidRad = fopen(fileRad);
Rad_File = fscanf(fidRad,'%g',[3 nbSteps]).';
fclose(fidRad);

bWantToSaveJPEG = true;
filename = strcat('radial',num2str(nbParticles),'.jpeg');

% r as a vector
r = [1 (nbSteps)];
for i = 1: (nbSteps)
    r(i) = Rad_File(i,1);
end

% Init of the vector holding the number of bound particles
nHist = [1 (nbSteps)];
for i = 1: (nbSteps)
    nHist(i) = Rad_File(i,2);
end

% And then we plot !
figure(3);
%hold on;
plot(r(:),nHist(:),'o','linestyle','-','color','b');
%hold off;
sTitle = strcat('Radial distribution of particles: N = ',num2str(nbParticles));
title(sTitle);
xlabel('r (ly)');
ylabel('N');
%axis([0, 0.16, -0.2, 1]);
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(3);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end    