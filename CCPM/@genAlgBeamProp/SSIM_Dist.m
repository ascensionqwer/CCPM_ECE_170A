function main
    % Main Script: Beam Propagation with Genetic Algorithm
    clear; clc; close all;

    %% Step 1: Define Parameters
    m = 1E3;          % Unit scaling (1 km in meters)
    um = 1E-6;        % Micrometer
    lambda = 1.55 * um; % Wavelength (1.55 micrometers)
    gridSize = 40;    % Computational grid size (arbitrary units)
    npts = 2^10;      % Number of grid points

    disp(['Initializing beam with parameters: ', ...
        'lambda = ', num2str(lambda), ', gridSize = ', num2str(gridSize), ...
        ', npts = ', num2str(npts)]);

    %% Step 2: Create Clean Gaussian Beam
    try
        % Use beamPropagation2D class
        cleanBeam = beamPropagation2D(lambda, gridSize, gridSize, npts, 'gauss');
        cleanBeamField = abs(cleanBeam.field_fList); % Intensity of the clean Gaussian beam
        save('cleanGaussianBeam.mat', 'cleanBeamField');
        disp('Clean Gaussian beam successfully created and saved.');
    catch ME
        disp(['Error creating clean Gaussian beam: ', ME.message]);
        return;
    end

    %% Step 3: Define Genetic Algorithm Parameters
    initVals = struct();
    initVals.herd_numHerd = 80;          % Population size
    initVals.herd_numGens = 100;         % Number of generations
    initVals.init_zProp = 100 * m;       % Propagation distance
    initVals.init_lambda = lambda;       % Wavelength
    initVals.init_gridSize = gridSize;   % Grid size
    initVals.const_npts = npts;          % Number of grid points
    initVals.init_beamType = 'gauss';    % Beam type

    disp('Genetic Algorithm parameters initialized.');

    %% Step 4: Initialize and Run Genetic Algorithm
    try
        % Initialize and run genAlgBeamProp class
        genAlg = genAlgBeamProp(initVals);
        genAlg.runAlg();
        disp('Genetic Algorithm successfully executed.');
    catch ME
        disp(['Error during Genetic Algorithm execution: ', ME.message]);
        return;
    end

    %% Step 5: Evaluate and Plot SSIM Performance
    try
        % Load clean Gaussian beam for reference
        load('cleanGaussianBeam.mat', 'cleanBeamField');
        cleanBeamFieldNorm = abs(cleanBeamField).^2 / max(abs(cleanBeamField(:)).^2);

        % Evaluate SSIM performance for each generation
        ssimHistory = NaN(initVals.herd_numGens, initVals.herd_numHerd);
        for gen = 1:initVals.herd_numGens
            disp(['Evaluating generation ', num2str(gen), '...']);
            for i = 1:length(genAlg.herd)
                if isfield(genAlg.herd(i), 'field_fList') && ~isempty(genAlg.herd(i).field_fList)
                    candidateField = abs(genAlg.herd(i).field_fList).^2;
                    candidateFieldNorm = candidateField / max(candidateField(:));
                    ssimHistory(gen, i) = SSIM_Dist(candidateFieldNorm, cleanBeamFieldNorm);
                    disp(['Generation ', num2str(gen), ', Candidate ', num2str(i), ': SSIM = ', num2str(ssimHistory(gen, i))]);
                else
                    disp(['Generation ', num2str(gen), ', Candidate ', num2str(i), ': field_fList is missing or empty.']);
                end
            end
        end

        % Plot SSIM Performance
        maxSSIMPerGeneration = max(ssimHistory, [], 2);
        figure;
        plot(1:initVals.herd_numGens, maxSSIMPerGeneration, '-o');
        title('SSIM Performance Across Generations');
        xlabel('Generation');
        ylabel('SSIM (Higher is Better)');
        grid on;

        disp('SSIM Performance successfully plotted.');
    catch ME
        disp(['Error during SSIM evaluation: ', ME.message]);
        return;
    end
end
