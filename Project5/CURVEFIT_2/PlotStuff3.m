fileToOpen = strcat('plotLeapfrog','.dat');
fidPosi = fopen(fileToOpen);
nbRows = 7;
Plot_File = fscanf(fidPosi,'%g',[20 nbRows]).';
fclose(fidPosi);
bWantToSaveJPEG = true;
filename = strcat('avgstd.jpeg');

% Init of our time vector
t = [1 (nbRows)];
for i = 1: (nbRows)
    t(i) = Plot_File(i,1);
end

% Init of the vector holding the number of bound particles
avg = [1 (nbRows)];
for i = 1: (nbRows)
    avg(i) = Plot_File(i,17);
end

std = [1 (nbRows)];
for i = 1: (nbRows)
    std(i) = Plot_File(i,18);
end

% And then we plot !
figure(1);
hold on;
plot(t(:),avg(:),'o','linestyle','-','color','k');
plot(t(:),std(:),'o','linestyle','-','color','r');
hold off;
title('Average distance to BOUND c.o.m. and std. dev.');
xlabel('N');
ylabel('average (black), std.dev. (red) (both in ly)');
axis([0, 2250, 0, 15]);
grid on;

% And eventually save the plot
if (bWantToSaveJPEG == true)  
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % To avoid 3D pictures
    imwrite(A,map,filename,'jpeg'); 
end