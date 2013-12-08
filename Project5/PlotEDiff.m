fileName = 'plotLeapfrog.dat';
fidPosi = fopen(fileName);
nbColumns = 19;
nbSimulation = 5;
Plot_File = fscanf(fidPosi,'%g',[nbColumns nbSimulation]).' ;
fclose(fidPosi);
bWantToSaveJPEG = true;
bFunctionOfEpsilon = false;
 
Plot_File(:,3);

% Init of our time vector
t = [1 nbSimulation];
for i = 1:nbSimulation
    t(i)= Plot_File((nbSimulation + 1- i),3);
end

t;
nbETotBefore =  [0 nbSimulation];
for i = 1 : nbSimulation
    nbETotBefore(i) = Plot_File((nbSimulation + 1- i),5);
end

nbETotBefore

% Init of our kinetic energy vector
nbETotAfter = [0 nbSimulation];
for i = 1: nbSimulation
    nbETotAfter(i) = Plot_File((nbSimulation + 1- i),6);
end
nbETotAfter

nbDiffEt = [0 nbSimulation];
for i = 1 : nbSimulation
    nbDiffEt(i) = nbETotBefore(i) - nbETotAfter(i);
end
nbDiffEt


% And then we plot !
%figure(1);
plot(t(:),nbDiffEt(:),'color',rand(1,3));
%title(['Plot for ' int2str(nbParticles) ' Particles']);
%xlabel('Time (in t_c_r_u_n_c_h)');
%ylabel('Number of bound particles');
grid on;
