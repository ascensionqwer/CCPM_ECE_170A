function main
    % Main Script: Compressed Results for Beam Propagation Experiments
    clear; clc; close all;

    %% Parameters
    m = 1e3;          % Base unit in meters
    cm = 1e-2 * m;
    um = 1e-6 * m;
    lambda = 1.55 * um; % Wavelength
    xSize = 4 * cm;     % Grid size
    ySize = 4 * cm;
    npts = 256;         % Grid points
    zDists = [50 * cm, 10 * m]; % Propagation distances

    %% Gaussian Beam Propagation
    disp('--- Gaussian Beam Propagation ---');
    beam_gauss = beamPropagation2D(lambda, xSize, ySize, npts, 'gauss');
    zR = pi * (xSize / 6)^2 / lambda; % Rayleigh range

    waistEstimates = [];
    waistTheory = [];
    for z = zDists
        propField = beam_gauss.forwardProp_FreeSpace2D(z);
        I = abs(propField).^2;

        % Estimate beam waist
        peakI = max(I(:));
        lineProfile = I(floor(npts / 2), :);
        w_est = sum(lineProfile > peakI * exp(-2)) * (xSize / npts) / 2; % Half max width
        w_theory = (xSize / 6) * sqrt(1 + (z / zR)^2);

        waistEstimates(end + 1) = w_est; %#ok<*AGROW>
        waistTheory(end + 1) = w_theory;

        disp(['z = ', num2str(z / cm), ' cm: Waist Estimate = ', num2str(w_est), ...
              ' mm, Theory = ', num2str(w_theory), ' mm']);
    end

    %% Laguerre-Gaussian Beam Propagation
    disp('--- Laguerre-Gaussian Beam Propagation ---');
    beam_lg = beamPropagation2D(lambda, xSize, ySize, npts, 'laguerre');
    propField_lg = beam_lg.forwardProp_FreeSpace2D(zDists(2));
    nmse = mean((abs(propField_lg(:)).^2 - abs(beam_lg.field_fList(:)).^2).^2) / ...
           mean(abs(beam_lg.field_fList(:)).^2);
    disp(['Laguerre-Gaussian NMSE: ', num2str(nmse)]);

    %% Sampling and Noise Robustness Studies
    disp('--- Sampling and Noise Robustness Studies ---');
    beam_small = beamPropagation2D(lambda, xSize, ySize, 64, 'gauss');
    propField_small = beam_small.forwardProp_FreeSpace2D(zDists(1));
    beam_large = beamPropagation2D(lambda, xSize, ySize, npts, 'gauss');
    propField_large = beam_large.forwardProp_FreeSpace2D(zDists(1));
    corrVal = corr2(abs(propField_small).^2, imresize(abs(propField_large).^2, [64, 64]));
    disp(['Grid Sampling Correlation: ', num2str(corrVal)]);

    % Add random noise to Gaussian beam
    beam_noise = beamPropagation2D(lambda, xSize, ySize, npts, 'gauss');
    beam_noise.field_fList = beam_noise.field_fList .* (1 + 0.1 * randn(size(beam_noise.field_fList)));
    propField_noise = beam_noise.forwardProp_FreeSpace2D(zDists(1));

    %% Plotting
    % Gaussian Beam Propagation at z = 50 cm and z = 10 m
    figure;
    subplot(1, 2, 1);
    imagesc(abs(beam_gauss.forwardProp_FreeSpace2D(zDists(1))).^2);
    title('Gaussian Beam (z = 50 cm)');
    xlabel('x (mm)'); ylabel('y (mm)');
    axis image; colormap jet; colorbar;
    c = colorbar;
    c.Label.String = 'Normalized Intensity';

    subplot(1, 2, 2);
    imagesc(abs(beam_gauss.forwardProp_FreeSpace2D(zDists(2))).^2);
    title('Gaussian Beam (z = 10 m)');
    xlabel('x (mm)'); ylabel('y (mm)');
    axis image; colormap jet; colorbar;
    c = colorbar;
    c.Label.String = 'Normalized Intensity';

    % Laguerre-Gaussian Beam at z = 10 m and Gaussian Beam with Noise
    figure;
    subplot(1, 2, 1);
    imagesc(abs(propField_lg).^2);
    title('Laguerre-Gaussian Beam (z = 10 m)');
    xlabel('x (mm)'); ylabel('y (mm)');
    axis image; colormap jet; colorbar;
    c = colorbar;
    c.Label.String = 'Normalized Intensity';

    subplot(1, 2, 2);
    imagesc(abs(propField_noise).^2);
    title('Gaussian Beam with Noise (z = 50 cm)');
    xlabel('x (mm)'); ylabel('y (mm)');
    axis image; colormap jet; colorbar;
    c = colorbar;
    c.Label.String = 'Normalized Intensity';

    % Sampling and Grid Comparison
    figure;
    subplot(1, 2, 1);
    imagesc(abs(propField_small).^2);
    title('Small Grid Sampling');
    xlabel('x (mm)'); ylabel('y (mm)');
    axis image; colormap jet; colorbar;
    c = colorbar;
    c.Label.String = 'Normalized Intensity';

    subplot(1, 2, 2);
    imagesc(imresize(abs(propField_large).^2, [64, 64]));
    title('Large Grid Resized to Match');
    xlabel('x (mm)'); ylabel('y (mm)');
    axis image; colormap jet; colorbar;
    c = colorbar;
    c.Label.String = 'Normalized Intensity';

    %% Display Results
    disp(['Waist Estimates: ', num2str(waistEstimates)]);
    disp(['Theoretical Waists: ', num2str(waistTheory)]);
    disp(['Laguerre-Gaussian NMSE: ', num2str(nmse)]);
    disp(['Grid Sampling Correlation: ', num2str(corrVal)]);
end
