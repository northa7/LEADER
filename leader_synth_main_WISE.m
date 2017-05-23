% -------------------------------------------------------------------------
%
% Recovering the shape elongation (p) and spin latitude (beta)
% distributions for an asteroid population using
% a brightness variation estimate (eta)
%
% NOTE: This routine utilizes synthetic data from DAMIT
% in order to compute the brightness data
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

p_tot = [];
beta_tot = [];
pb_tot = [];

% DAMIT file listing
fid = fopen('asteroideja.txt');
temp = textscan(fid, '%s');
fclose(fid);
damit = temp{1};

% Choose a single peak for the (p, beta) distribution
% (random, unless already declared)
if ( ~exist('P_PEAK','var') || ~exist('B_PEAK', 'var') )
    P_PEAK = 0.6*rand+0.35;
    B_PEAK = 1.5*rand+0.05;
end
% Fix P_PEAK and B_PEAK, if needed
%P_PEAK = 0.7;
%B_PEAK = 0.2;

% The number of asteroids to be analysed
if ( ~exist('Nkierroksia','var') || isempty(Nkierroksia) )
    Nkierroksia = 1000;
end
Napparlist = zeros(Nkierroksia, 1);

% -------------
% Forward model
% -------------

for kierroksia=1:Nkierroksia
    p_test = 0;
    while ( abs(p_test-P_PEAK) > 0.075 )
        tiedostonimi = damit{ ceil(length(damit)*rand) };

        [X,Y,Z,F] = damit_model(tiedostonimi, false);    
        R = [X,Y,Z];
        % Stretching (random, with stretching factor being at least 1)
        stretch = diag(max([1; 1; 1], 2*abs(randn(3,1))));
        R = R*stretch*eye(3);

        % Go through the faces to obtain the elongation p = b/a
        leader_ellipsoid
        p_test = p;
        
        % With a 10% chance, accept a random p (as long as p > 0.45)
        if ( rand > 0.9 && p_test > 0.45 )
            break
        end
    end
    
    % Display progress (disable, if the feature is not wanted)
    if ( mod(kierroksia,50) == 0)
        disp([num2str(kierroksia)])
    end
    
    % Read the data files to acquire brightness data, and compute etas
    if ( exist('lcg_files', 'var') && Nkierroksia == length(lcg_files) )
        % Whole database => go through the asteroids in order
        % (otherwise the order is random)
        mones = kierroksia;
    end
    leader_brightness_synth_WISE

    % Keep information on the number of apparitions
    Napparlist(kierroksia) = Nappar;
    
    % Add to the list of amplitudes
    A_tot = [A_tot; A];
    p_tot = [p_tot; p];
    beta_tot = [beta_tot; beta];
    pb_tot = [pb_tot; p, beta];
end

% Construct the CDF of A
Asort = sort(A_tot);
CDFA = linspace(1/length(Asort), 1, length(Asort))';

p = p_tot;
beta = beta_tot;

% Draw the marginal DFs of the actual p and beta
if ( length(p) > 1 )
    weightbar(p, ones(size(p)))
    xlabel('p'), ylabel('n')
    drawnow
end
if ( length(beta) > 1 )
    weightbar(beta, ones(size(beta)), 0, pi/2, 16)
    xlabel('\beta'), ylabel('n')
    drawnow
end
if ( length(p) >= 1 && length(beta) >= 1 )
    % Contour plot of the forward model
    figure
    [w, cc] = hist3(pb_tot, ...
        {linspace(0.05, 0.95, 10), linspace(0.05, 1.55, 16)});
    [PTOT,BTOT]=meshgrid(cc{1}, cc{2});
    contourf(PTOT,BTOT,w')
    colorbar
    xlabel('p'), ylabel('\beta')
    drawnow
end

% Solve the DFs of p and beta from the inverse problem
leader_invert

% Solution plots
leader_plots

% Visual postprocessing for the solution contour
leader_postprocess_WISE