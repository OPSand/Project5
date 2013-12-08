fileLeap = strcat('plotLeapfrog.dat');
fidLeap = fopen(fileLeap);
Plot_File = fscanf(fidLeap,'%g',[18 4]).';
fclose(fidLeap);

simID = 3;
nbParticles = 2000;
nbSteps = 67;
fileRad = strcat('radial_after_',simID,'_leapfrog.dat');
fidRad = fopen(fileRad);
Rad_File = fscanf(fidRad,'%g',[3, nbSteps]).';
fclose(fidRad);

bWantToSaveJPEG = true;
filename = strcat('fit',nbParticles,'.jpeg');

% r as a vector
r = [1 (nbSteps)];
for i = 1: (nbSteps)
    r(i) = Plot_Rad(i,1);
end

n0 = Plot_File(4,11);
r0 = Plot_File(4,12);

% Init of the vector holding the number of bound particles
n = [1 (nbSteps)];
for i = 1: (nbSteps)
    n(i) = (n0 / (1 + (r(i)/r0)^4));
end

% Init of the vector holding the number of bound particles
nbBounds2 = [1 (nbRows)];
for i = 1: (nbRows)
    nbBounds2(i) = Plot_File(i,8);
end

% And then we plot !
figure(1);
hold on;
plot(r(:),nbBounds(:),'o','linestyle','-','color','r');
plot(r(:),nbBounds2(:),'o','linestyle','-','color','b');
hold off;
title(['Energy Conservation']);
xlabel('Epsilon (ly)');
ylabel('Relative change in energy (1 = 100%)');
axis([0, 0.16, -0.2, 1]);
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end
    