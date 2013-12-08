fileToOpen = strcat('plotLeapfrog','.dat');
fidPosi = fopen(fileToOpen);
nbRows = 16;
Plot_File = fscanf(fidPosi,'%g',[18 nbRows]).';
fclose(fidPosi);
bWantToSaveJPEG = true;
filename = strcat('boundFraction.jpeg');

% Init of our time vector
t = [1 (nbRows)];
for i = 1: (nbRows)
    t(i) = Plot_File(i,2);
end

% Init of the vector holding the number of bound particles
nbBounds = [1 (nbRows)];
for i = 1: (nbRows)
    nbBounds(i) = Plot_File(i,18);
end

% And then we plot !
figure(1);
%hold on;
plot(t(:),nbBounds(:),'o','linestyle','-','color','k');
%hold off;
title(['Bound particles']);
xlabel('Epsilon (ly)');
ylabel('Fraction of bound particles');
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end
    