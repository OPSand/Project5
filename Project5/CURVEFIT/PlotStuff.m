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
filename = strcat('fit',num2str(nbParticles),'.jpeg');
filename2 = strcat('fit',num2str(nbParticles),'_2.jpeg');

% r as a vector
r = [1 (nbSteps)];
rN13 = [1 (nbSteps)];
for i = 1: (nbSteps)
    r(i) = Rad_File(i,1);
    rN13(i) = r(i) / (nbParticles^(-1/3));
end

n0 = Plot_File(4,12);
r0 = Plot_File(4,13);

n = [1 (nbSteps)];
nN2 = [1 (nbSteps)];
for i = 1: (nbSteps)
    n(i) = (n0 / (1 + (r(i)/r0)^4));
    nN2(i) = n(i) / (nbParticles^2);
end

% Init of the vector holding the number of bound particles
nHist = [1 (nbSteps)];
nHistN2 = [1 (nbSteps)];
for i = 1: (nbSteps)
    nHist(i) = Rad_File(i,3);
    nHistN2(i) = nHist(i) / (nbParticles^2);
end

% And then we plot !
figure(3);
semilogy(r(:),n(:),'o','linestyle','-','color','r')
hold on;
semilogy(r(:),nHist(:),'o','linestyle','-','color','b');
hold off;
sTitle = strcat('Curve fitting: n0 = ',num2str(n0),', r0 = ',num2str(r0), ', ',num2str(nbParticles),' particles');
title(sTitle);
xlabel('r (ly)');
ylabel('n (particles/ly³');
%axis([0, 0.16, -0.2, 1]);
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(3);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end

% And then we plot !
figure(4);
semilogy(rN13(:),nN2(:),'o','linestyle','-','color','r')
hold on;
semilogy(rN13(:),nHistN2(:),'o','linestyle','-','color','b');
hold off;
sTitle = strcat('Curve fitting: n0 = ',num2str(n0),', r0 = ',num2str(r0), ', ',num2str(nbParticles),' particles');
title(sTitle);
xlabel('r/N^-^1^/^3 (ly)');
ylabel('n/N^2 (1/ly³');
%axis([0, 0.16, -0.2, 1]);
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(4);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end
    