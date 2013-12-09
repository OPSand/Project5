fileName = 'plotLeapfrog.dat';
fidPosi = fopen(fileName);
nbColumns = 19;
nbSimulation = 4;
Plot_File = fscanf(fidPosi,'%g',[nbColumns nbSimulation]).' ;
fclose(fidPosi);
bWantToSaveJPEG = true;

Plot_File

% Init of our time vector
t = [1 nbSimulation];
for i = 1:nbSimulation
    t(i)= Plot_File((nbSimulation + 1- i),3);
end
t;

t; % Remove the ";" to know what's inside me !

nbETotBefore =  [0 nbSimulation];
for i = 1 : nbSimulation
    nbETotBefore(i) = Plot_File((nbSimulation + 1- i),5);
end

nbETotBefore; % Remove the ";" to know what's inside me !

% Init of our kinetic energy vector
nbETotAfter = [0 nbSimulation];
for i = 1: nbSimulation
    nbETotAfter(i) = Plot_File((nbSimulation + 1- i),6);
end
nbETotAfter; % Remove the ";" to know what's inside me !

nbDiffEt = [0 nbSimulation];
for i = 1 : nbSimulation
    nbDiffEt(i) = nbETotBefore(i) - nbETotAfter(i);
end

nbDiffEt; % Remove the ";" to know what's inside me !


% And then we plot !
%figure(1);
% plot(t(:),nbDiffEt(:),'o','lineStyle','-','color',rand(1,3)) %color',rand(1,3),
% title('Plot for Leapfrog -- Conservation of Energy -- Epsilon = 0,1');
% xlabel('Size time step (in t_c_r_u_n_c_h)');
% ylabel('E_i_n_i_t - E_f_i_n_a_l');
% grid on;
% 
% % And eventually save the plot
% if (bWantToSaveJPEG == true)  
%     saveName = 'Plot_Conservation_E_Leapfrog.jpeg';
%     frame = getframe(1);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256); % To avoid 3D pictures
%     imwrite(A,map,saveName,'jpeg'); 
% end

