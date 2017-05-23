% -------------------------------------------------------------------------
%
% Recovering the shape elongation (p) and spin latitude (beta)
% distributions for an asteroid population using
% a brightness variation estimate (eta)
%
% - Shape elongation: 0 < p <= 1, where 0 means elongated, 1 means
% spherical; p=b/a (for an ellipsoid with semiaxes a >= b >= c)
% - Spin latitude: 0 <= beta <= pi/2, where beta is measured from the pole
% (0: perpendicular to the ecliptic plane, pi/2: in the ecliptic plane)
%
% NOTE: You must have preloaded an Nx1-sized cell variable
% called 'lcg_files' which contains the directory paths of
% all the data files. If the variable isn't found, the
% routine attempts to load a fixed .mat file that contains
% an example of the directory listing variable 'lcg_files'.
%
% This version uses WISE as the database
%
% -------------------------------------------------------------------------

A_tot = [];
Npisteita = [];
marginplots = true;

% The number of asteroids to be analysed
if ( ~exist('Nkierroksia','var') || isempty(Nkierroksia) )
    Nkierroksia = 1000;
end
Napparlist = zeros(Nkierroksia, 1);

for kierroksia=1:Nkierroksia
    % Display progress (disable, if the feature is not wanted)
    if ( mod(kierroksia,1000) == 0)
        disp([num2str(kierroksia)])
    end
    
    % Read the data files to acquire brightness data, and compute etas
    if ( exist('lcg_files', 'var') && Nkierroksia == length(lcg_files) )
        % Whole database => go through the asteroids in order
        % (otherwise the order is random)
        mones = kierroksia;
    end
    lcg_read_WISE
    
    % Keep information on the number of apparitions
    Napparlist(kierroksia) = Nappar;

    % Add to the list of amplitudes
    A_tot = [A_tot; A];
end

% Construct the CDF of A
Asort = sort(A_tot);
CDFA = linspace(1/length(Asort), 1, length(Asort))';

% Solve the DFs of p and beta from the inverse problem
leader_invert

% Solution plots
leader_plots

% Visual postprocessing for the solution contour
leader_postprocess_WISE