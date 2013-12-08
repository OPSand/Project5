fileToOpen = strcat('plotLeapfrog','.dat');
fidPosi = fopen(fileToOpen);
nbRows = 3;
Plot_File = fscanf(fidPosi,'%g',[20 nbRows]).';
fclose(fidPosi);
bWantToSaveJPEG = true;
filename = strcat('fitscaling_n0.jpeg');
filename2 = strcat('fitscaling_r0.jpeg');

% Init of our time vector
t = [1 (nbRows)];
for i = 1: (nbRows)
    t(i) = Plot_File(i,1);
end

% Init of the vector holding the number of bound particles
nbBounds = [1 (nbRows)];
for i = 1: (nbRows)
    nbBounds(i) = Plot_File(i,15);
end

rbBounds = [1 (nbRows)];
for i = 1: (nbRows)
    rbBounds(i) = Plot_File(i,16);
end

% And then we plot !
figure(1);
hold on;
plot(t(:),nbBounds(:),'o','linestyle','-','color','k');
hold off;
title('Scaling of n0 with N');
xlabel('N');
ylabel('n0 / N²');
%axis([0, 2250, 0, 0.01]);
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end

% And then we plot !
figure(2);
hold on;
plot(t(:),rbBounds(:),'o','linestyle','-','color','k');
hold off;
title('Scaling of r0 with N');
xlabel('N');
ylabel('r0 / N^(^-^1^/^3^)');
%axis([0, 2250, 0, 1]);
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(2);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename2,'jpeg'); 
end