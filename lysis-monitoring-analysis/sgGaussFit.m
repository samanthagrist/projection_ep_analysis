%**************************************************************************%
%*                                                                        *%
%*                             Gaussian Fit                               *%
%*                                                                        *%
%*                           Samantha M. Grist                            *%
%*                           v0.1 October 2018                            *%
%*                                                                        *%
%**************************************************************************%

function [a, mu, sigma, rsq] = sgGaussFit(x, y)

%Performs a Gaussian fit to x-y data
%INPUTS
% x: numerical array of x data values to fit
% y: numerical array of x data values to fit

% OUTPUTS
% a: Gaussian peak amplitude
% mu: Gaussian peak centre
% sigma: Gaussian peak width sigma
% rsq: R-squared value for the fit

% fit each ROI cross-section to a Gaussian to find the peak location
fit_type = 'gauss1';
fit_options = fitoptions(fit_type); 


% Set the sigma bounds
sigma_min = 0;
%sigma_max = 2*length(y); SG changed for sim data 2018-11-22
sigma_max = 10*max(x);
% estimate peak width
[maxVal, idx] = max(y);
for i = idx:length(y)
    if y(i) < 0.6*max(y)
        rightBound = i;
        break;
    end
end
for i = 1:idx
    if y(i) > 0.6*max(y)
        leftBound = i;
        break;
    end
end
sigma_start = (x(rightBound)-x(leftBound))/2;

% Set the peak center bounds
x_min = min(x);
x_max = max(x);

% Set the ampitude bounds
a_min = 0;
a_max = 5*max(y);

% set the upper and lower bounds. correct for difference in c and
% sigma terms
fit_options.Lower = [a_min, x_min, (sigma_min * sqrt(2))];
fit_options.Upper = [a_max, x_max, (sigma_max * sqrt(2))];
fit_options.StartPoint = [max(y), 0, (sigma_start * sqrt(2))];
fit_options.Robust = 'LAR';

% Fit the peaks
[fit_object, gof] = fit(x, y, fit_type, fit_options);
figure(300); plot(fit_object, x, y);

% Get the coefficients
fit_coeffs = coeffvalues(fit_object);
rsq = gof.rsquare;
a = fit_coeffs(1);
mu = fit_coeffs(2);
sigma = fit_coeffs(3)/sqrt(2);
                    
             