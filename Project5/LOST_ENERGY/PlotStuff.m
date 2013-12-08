fileLeap = strcat('plotLeapfrog.dat');
fidLeap = fopen(fileLeap);
Plot_File = fscanf(fidLeap,'%g',[20 4]).';
fclose(fidLeap);

bWantToSaveJPEG = true;
filename = strcat('lostenergy.jpeg');

n = [1 (nbSteps)];
for i = 1: (nbSteps)
    n(i) = Plot_File(i, 1);
end

dE = [1 (nbSteps)];
for i = 1: (nbSteps)
    dE(i) = Plot_File(i, 20);
end

% And then we plot !
figure(1);
plot(n(:),dE(:),'o','linestyle','-','color','k')
sTitle = strcat('Effect of N on energy lost due to particle ejection');
title(sTitle);
xlabel('N');
ylabel('Energy lost (solar masses*ly²/tcrunch²)');
%axis([0, 0.16, -0.2, 1]);
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end
    