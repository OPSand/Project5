fileName = 'plotLeapfrog.dat';
fidPosi = fopen(fileName);
nbColumns = 18;
nbSimulation = 5;
Plot_File = fscanf(fidPosi,'%g',[nbColumns nbSimulation]).' ;
fclose(fidPosi);
bWantToSaveJPEG = true;
bFunctionOfEpsilon = false;
 

% Init of our kinetic energy vector
nbBoundsEk = [0 nbSimulation];
for i = 1: nbSimulation
    nbBoundsEk(i) = Plot_File(i,10);
end

nbBoundsEk; % Decomment to check that everything is well behaving

% Init of our kinetic energy vector
 nbBoundsEp = [0 nbSimulation];
 for i = 1: nbSimulation
     nbBoundsEp(i) = Plot_File(i,11);
 end

 nbBoundsEp; % Decomment to check that everything is well behaving

nbBoundsVirial = [0 nbSimulation];
for i = 1: nbSimulation
    nbBoundsVirial(i) = nbBoundsEp(i)/nbBoundsEk(i);
end

nbBoundsVirial; % Decomment to check that everything is well behaving


if bFunctionOfEpsilon == true

    % Init of our time vector
    t = [0 (nbSimulation)];
    for i = 1:nbSimulation
        t(i)= Plot_File(i,2);
    end

    
    % And then we plot !
    figure(1);
    plot(t(:),nbBoundsVirial(:),'o','lineStyle','-');
    title('Ration E_p/E_k in function of epsilon');
    xlabel('Epsilon(ly)');
    ylabel('E_p/E_k');
    grid on;

    % And eventually save the plot
    if (bWantToSaveJPEG == true)  
        saveName = 'Plot_Virial_Epsilon.jpeg';
        frame = getframe(1);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256); % To avoid 3D pictures
        imwrite(A,map,saveName,'jpeg'); 
    end
% We want if in function of the number of particles!     
else
    nbMaxParticles = 300;
    nbMinParticles = 100;
    % Init of our nbParticles vectors
    % On va de 100 à 300 en 6 étapes
    nbParticles = [0 nbSimulation];
    for i = 1: nbSimulation
        nbParticles(i) = nbMinParticles + (i-1)*(nbMaxParticles - nbMinParticles)/(nbSimulation -1);
    end
    nbParticles;  % Decomment to check that everything is well behaving
    
    % And then we plot !
    figure(1);
    plot(nbParticles(:),nbBoundsVirial(:),'o','lineStyle','-');%,'Color','g','+');
    title('Ration E_p/E_k in function of the number of particles');
    xlabel('N');
    ylabel('E_p/E_k');
    grid on;

    % And eventually save the plot
    if (bWantToSaveJPEG == true)  
        saveName = 'Plot_Virial_N.jpeg';
        frame = getframe(1);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256); % To avoid 3D pictures
        imwrite(A,map,saveName,'jpeg'); 
    end
    
end
   