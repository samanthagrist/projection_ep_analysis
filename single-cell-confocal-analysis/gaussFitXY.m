%*************************************************************************%
%**                    Gaussian fit XY peak widths                      **%
%*************************************************************************%
% Samantha M. Grist
% v.0.2 - 2019-07-08
% for simple Z profile analysis
 
function [fit_object, gof] = gaussFitXY(x,y)
    fit_type = 'gauss1';
    fit_options = fitoptions(fit_type); 

    % Set the sigma bounds
    sigma_min = 0;
    sigma_max = 2*length(y);

    % Set the peak center bounds
    x_min = min(x);
    x_max = max(x);

    % Set the ampitude bounds
    a_min = 0;
    if max(y) > 0
        a_max = 1.5*max(y);
    else % if the whole range is below zero, artificially shift it up so we can fit.
        y = y - min(y);
        a_max = 1.5*max(y);
    end

    % set the upper and lower bounds. correct for difference in c and
    % sigma terms
    fit_options.Lower = [a_min, x_min, (sigma_min * sqrt(2))];
    fit_options.Upper = [a_max, x_max, (sigma_max * sqrt(2))];
%                 fit_options.Lower
%                 fit_options.Upper

    % Fit the peaks
    [fit_object, gof] = fit(x, y, fit_type, fit_options);
end