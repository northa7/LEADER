% Read WISE data files (code meant for synthesized version)
if ( ~exist('lcg_files', 'var') )
    load lcg_files_WISE
end

% If the object we want to read hasn't been declared separately,
% choose a random object from the database
if ( ~exist('mones', 'var') || isempty(mones) )
    mones = ceil(length(lcg_files)*rand);
end
tiednimi = lcg_files{mones};
fid = fopen(tiednimi);
clear mones

% The first row has the number observations
temp = textscan(fid, '%f', 1, 'Delimiter', '\n');
norows = temp{1};
dates = zeros(norows, 1);
e_sun = zeros(norows, 3);
e_earth = zeros(norows, 3);
flux = cell(norows, 4);
fluxerr = cell(norows, 4);
for i=1:norows
    % First, date and the number of filters
    temp = textscan(fid, '%f %f', 1, 'Delimiter', '\n');
    % Error check: break, if norows > data rows
    if ( isempty(temp{1}) )
        break
    end
    dates(i) = temp{1};
    nofilters = temp{2};
    
    % Next, the direction of the Sun
    temp = textscan(fid, '%f %f %f', 1, 'Delimiter', '\n');
    e_sun(i,:) = [temp{1}, temp{2}, temp{3}];
    
    % Followed by the direction of the Earth
    temp = textscan(fid, '%f %f %f', 1, 'Delimiter', '\n');
    e_earth(i,:) = [temp{1}, temp{2}, temp{3}];
    
    % Filters
    for j=1:nofilters
        temp = textscan(fid, '%f %f %f %f', 1, 'Delimiter', '\n');
        % Filter code (0-3)
        filter = temp{4};
        % The first number is the wavelength
        % The second and third numbers are flux and the error of flux
        flux{i, filter+1} = temp{2};
        fluxerr{i, filter+1} = temp{3};
    end
    
    % Skip 2 blank rows
    temp = textscan(fid, '%c', 1, 'Delimiter', '\n');
end
fclose(fid);

% Construct vectors from flux and flux error
flux1 = [];
flux2 = [];
flux3 = [];
flux4 = [];
flux1e = [];
flux2e = [];
flux3e = [];
flux4e = [];

for i=1:norows
    if ( ~isempty(flux{i,1}) )
        flux1 = [flux1; flux{i,1}];
        flux1e = [flux1e; fluxerr{i,1}];
    end
end

for i=1:norows
    if ( ~isempty(flux{i,2}) )
        flux2 = [flux2; flux{i,2}];
        flux2e = [flux2e; fluxerr{i,2}];
    end
end

for i=1:norows
    if ( ~isempty(flux{i,3}) )
        flux3 = [flux3; flux{i,3}];
        flux3e = [flux3e; fluxerr{i,3}];
    end
end

for i=1:norows
    if ( ~isempty(flux{i,4}) )
        flux4 = [flux4; flux{i,4}];
        flux4e = [flux4e; fluxerr{i,4}];
    end
end

% Compute total errors
flux_tot_err = [norm(flux1e), norm(flux2e), norm(flux3e), norm(flux4e)];
% Discard empty filters
flux_tot_err(flux_tot_err == 0) = Inf;
% Best filter
bestf = find(flux_tot_err == min(flux_tot_err), 1);

% Collect the intensity data from the best filter
L_big = [];
indeksit = [];
for i=1:norows
    if ( ~isempty(flux{i,bestf}) )
        L_big = [L_big; flux{i, bestf}];
        indeksit = [indeksit; i];
    end
end
% Take only the observations done with this filter
e_sun = e_sun(indeksit, :);
e_eath = e_earth(indeksit, :);
dates = dates(indeksit);

% Normalize
e_sun = normr(e_sun);
e_earth = normr(e_earth);

% Remove the measurements where the phase angle > ang_tol
ang_tol = deg2rad(30);
poistettavat = false(length(L_big),1);
ang = zeros(length(L_big), 1);
for i=1:length(L_big)
    ang(i) = acos( e_sun(i,:)*e_earth(i,:)' );
    if ( ang(i) > ang_tol )
        poistettavat(i) = true;
    end
end
dates(poistettavat) = [];
L_big(poistettavat) = [];
e_sun(poistettavat, :) = [];
e_earth(poistettavat, :) = [];
ang(poistettavat) = [];