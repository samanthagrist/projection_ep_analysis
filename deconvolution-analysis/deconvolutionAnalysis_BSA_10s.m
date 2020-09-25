%*************************************************************************%
%*************************************************************************%
%**                                                                     **%
%**                                                                     **%
%**                  DIFFUSIONAL PSF DECONVOLUTION                      **%
%**                       October 11th, 2019                            **%
%**                        Samantha M. Grist                            **%
%**                                                                     **%
%**                                                                     **%
%*************************************************************************%
%*************************************************************************%

% Assess the quality of physics-driven postprocessing to reduce the
% x-y peak width of separated protein bands

close all;
clear;
clc;

% set up the user configurable parameters
protein = {'BSA'};
MW = 66.5; %[66.5, 42.7, 150];
Rh = 3.011185362; %[3.011185362, 2.669790732, 3.824816335]; % from diffusion calculation spreadsheet
gels = 7;
gelText = {'7%T'};
Tstart = 4; % starting temperature, C
filename = '7T_10s-2_Z37_BSA_80um.tif';
t = 20; % diffusion time 7T_5s-1 t=15; 7T_10s-2 t=20s; 7T_15s-2 t=26s
filenameOut = '7T_10s-2_Z37_BSA_80um_Processed';
xi = 16e-6; % initial band size (m)
imgScale = 0.339; % microns/pixel
xsecWidth = 100; % width of the cross section to sum and analyze (microns)
numXsecs = 3;

% set up the auto parameters
Tstart = Tstart+273.15;
kB = (1.38*10^-23);
colours = 'kbgrcm';
linestylelist = [{':'}, {'--'}, {'-.'}, {''}];
C = 1e-6; % initial protein concentration in injected band
% estimate that the diffusion time also includes ~half the photocapture time
%t = t+22;
kT = kB*Tstart;
mu = ((2.761*10^-3)*exp(1713/Tstart))/1000;



% set up colormap
cyanLut = [(1:(-1/255):0)',(1:(-0.2/255):0.8)',(1:(-0.2/255):0.8)']; % 0 is white; max is higher vis cyan
magentaLut = [(1:(-0.2/255):0.8)',(1:(-1/255):0)',(1:(-0.2/255):0.8)']; % 0 is white; max is higher vis magenta
redLut = [(1:(-0.2/255):0.8)',(1:(-1/255):0)',(1:(-1/255):0)']; % 0 is white; max is higher vis red
imgLUT = cyanLut;


% set up the point spread function
% set up the x, y, z parameters (in m except sampling)
sampling = imgScale; % um
x = (-100:sampling:100).*1e-6;
y = (-100:sampling:100).*1e-6;
z = (-100:sampling:100).*1e-6;
[X, Y, Z,] = meshgrid(x, y, z);
[~, xindex] = min(abs(x));
[~, yindex] = min(abs(y));
[~, index] = min(abs(z));

% find the diffusion coefficients
Dsol = (kT/(6*3.14*mu*Rh*10^-9));
gelSolnRatio = exp(-3.03*((Rh*10)^0.59)*(gels/100)^0.94);
Dgel = Dsol*gelSolnRatio;

% set up the psf for the diffusion time of interest
if t == 0 % at t=0 just have a point
    psf_diff = zeros(length(x), length(y), length(z));
    psf_diff(xindex, yindex, index) = 1; % initial concentration is 1/size of voxel
else % t>0 use the diffusion function.  % initial injected mass is 1[mol/L]*size of voxel
    psf_diff(:, :, :) = 1*(sampling*1e-6)^3./((4*pi().*t.*Dgel).^(3/2)).*exp(-(X.^2+Y.^2+Z.^2)./(4.*Dgel.*t));
end

save([filenameOut,'.mat']);

 % find the psf x-y profile at the z-centre of the psf
temp = squeeze(psf_diff(:, :, index));
% normalize the psf so it sums to 1
psfslice = (double(temp)/double(sum(temp(:))));

% read in the image
I = imread(filename);

%% 

       
I_disp = sgautoscale(I);
hOrig = figure; imshow(I_disp, imgLUT); title('Original image');
figure; imshow(sgautoscale(psfslice)); title('PSF');

% run the deconvolution
J = deconvlucy(I, psfslice);
figure; imshow(sgautoscale(J), imgLUT); title('Deconvolved image');

% save images with new colormap, with same display parameters
% deconvolved image should have higher intensity, so scale both to its max
J_max = max(J(:));
Jscaled = sgautoscale(J);
Iscaled = uint8(double(255*(double(I)./double(J_max))));
imwrite(J, imgLUT, [filenameOut, '_deconv.tif']);
imwrite(I, imgLUT, [filenameOut, '_scaled.tif']);
imwrite(J, imgLUT, [filenameOut, '_deconv.png'], 'BitDepth', 8);
imwrite(I, imgLUT, [filenameOut, '_scaled.png'], 'BitDepth', 8);

save([filenameOut,'.mat']);

%% 

% Make the intensity profiles and find the peaks
ADisp = I_disp;
for i=1:numXsecs
    % Allow the user to pick the 3 cross-sections to sum for Gaussian fitting
    figure(hOrig);
    title('Click in the centre of the ROI (horizontal rectangle) that you wish to analyze')
    [x,y,button] = ginput(1);
    % Burn the rectangle into the display image
    xSecWPix = xsecWidth/imgScale;
    ADisp = insertShape(ADisp, 'Rectangle', [1,y-xSecWPix/2,size(ADisp, 1)-1,xSecWPix], 'Color', 'r');
    pos = get(hOrig, 'Position');
    figure(hOrig); imshow(ADisp), 
    set(hOrig, 'Position', pos);
    
    % make the cross section
    ASub = I(y-xSecWPix/2:y+xSecWPix/2, 1:size(ADisp, 1));
    JSub = J(y-xSecWPix/2:y+xSecWPix/2, 1:size(ADisp, 1));
    figure(300); imshow(sgautoscale(ASub));
    xsec(i, :) = sum(ASub, 1)';
    xsecDeconv(i, :) = sum(JSub, 1)';
    
    % find the peaks (denoise with Gaussian filter first to minimize extra peaks)
    w = gausswin(50);
    [pks, locs] = findpeaks(filter(w, 1, xsec(i, :)), 'SortStr', 'descend', 'NPeaks', 3);
    
    % remove any peaks close to the edges of the image
    rowsToDelete = (locs < xsecWidth) | (locs > size(ADisp, 1)-xsecWidth);
    locs(rowsToDelete) = [];
    pkLocs(i) = {locs};
    
    % Gaussian fit the peaks
    for pkIdx = 1:length(locs)
        disp(['Peak ', num2str(pkIdx)]);
        pkLoc = locs(pkIdx);
        x = ((pkLoc-xSecWPix/2):(pkLoc+xSecWPix/2))'.*imgScale; % x in microns
        y = xsec(i, uint16(pkLoc-xSecWPix/2):uint16(pkLoc+xSecWPix/2))';
        yDeconv = xsecDeconv(i, uint16(pkLoc-xSecWPix/2):uint16(pkLoc+xSecWPix/2))';
        
        fit_type = 'gauss1';
        fit_options = fitoptions(fit_type); 
        fit_optionsD = fitoptions(fit_type); 

        % Set the sigma bounds
        sigma_min = 0;
        sigma_max = 2*length(y);

        % Set the peak center bounds
        x_min = min(x);
        x_max = max(x);

        % Set the ampitude bounds
        a_min = 0;
        a_max = 1.5*max(y);
        a_maxD = 1.5*max(yDeconv);

        % set the upper and lower bounds. correct for difference in c and
        % sigma terms
        fit_options.Lower = [a_min, x_min, (sigma_min * sqrt(2))];
        fit_options.Upper = [a_max, x_max, (sigma_max * sqrt(2))];
        fit_optionsD.Lower = [a_min, x_min, (sigma_min * sqrt(2))];
        fit_optionsD.Upper = [a_maxD, x_max, (sigma_max * sqrt(2))];

        % Fit the peaks
        [fit_object, gof] = fit(x, y, fit_type, fit_options);
        [fit_objectD, gofD] = fit(x, yDeconv, fit_type, fit_optionsD);
        figure(300); plot(fit_object, x, y); 
        figure(301); plot(fit_objectD, x, yDeconv); pause(0.01);

        % save the amplitude, mean, peak width, and goodness of fit
        fit_coeffs = coeffvalues(fit_object);
        fit_coeffsD = coeffvalues(fit_objectD);
        % amplitude
        a(i, pkIdx) = fit_coeffs(1);
        aD(i, pkIdx) = fit_coeffsD(1);
        % mu
        mu(i, pkIdx) = fit_coeffs(2);
        muD(i, pkIdx) = fit_coeffsD(2);
        % sigma
        sigma(i, pkIdx) = fit_coeffs(3)/sqrt(2);
        sigmaD(i, pkIdx) = fit_coeffsD(3)/sqrt(2);
        % rsq
        rsq(i, pkIdx) = gof.rsquare;
        rsqD(i, pkIdx) = gofD.rsquare;
        disp(['R2 = ', num2str(rsq(i, pkIdx)), ', R2 Deconv = ', num2str(rsqD(i, pkIdx))]);
        
        % find the AUCs for each peak
        [~,pkCentre] = min(abs(x-mu(i, pkIdx)));
        [~,locIdx] = min(abs(x-double(pkLoc*imgScale)));
        absPkCentre = uint16(pkLoc+(pkCentre-locIdx));
        idxSigma = uint16(sigma(i, pkIdx)/imgScale);
        pkToSum = xsec(i, (absPkCentre-2*idxSigma):(absPkCentre+2*idxSigma));
        pkToSumD = xsecDeconv(i, (absPkCentre-2*idxSigma):(absPkCentre+2*idxSigma));
        AUC(i, pkIdx) = sum(pkToSum);
        AUCD(i, pkIdx) = sum(pkToSumD);
    end
end

save([filenameOut,'.mat']);

deltaMu = (mu-muD);
deltaSigma = (sigma-sigmaD);
deltaAUC = (AUC-AUCD);



%% 
for i=1:numXsecs
        % plot the cross sections
    position = [1:length(xsec(i, :))]*imgScale;
    h = figure(88);
    hold on;
    plot(position, xsec(i,:), 'Color', [0, 0.8, 0.8]);
    plot(position, xsecDeconv(i,:), 'Color', [0, 0.4, 0.4]);
    xlabel('Position (\mum)');
    ylabel('AFU');
end

saveas(h, [filenameOut, '_profiles.fig']);

