% Read WISE data files

% Load the WISE database
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

% The first row has the number of observations
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

% For each eta, the data must be measured
% within a time span set by date_tol
date_tol = 60;
% The number of points wanted for an eta estimate
wanted = 5;
A = [];
i = 1;
% Backup the data
dates_back = [];
ang_back = [];
temp_kulma = [];
L_back = [];
pointsperapp = [];
Nappar = 0;
while ( i < length(L_big) )
    L = L_big(i);
    for j=i+1:length(L_big)
        % Accept the data if the time span <= date_tol, quit otherwise
        if ( dates(j) - dates(i) <= date_tol )
            L = [L; L_big(j)];
            % Close the loop in case it's the last measurement of the LC
            if ( j == length(L_big) )
                i_old = i;
                i = j;
            end
        else
            i_old = i;
            i = j;
            break
        end
    end
    
    % Store the data
    dates_back = [dates_back; dates(i_old:i_old+length(L)-1)];
    ang_back = [ang_back; ang(i_old:i_old+length(L)-1)];
    L_back = [L_back; L];
    pointsperapp = [pointsperapp; length(L)];
    Nappar = Nappar + 1;
    
    temp = ang(i_old:i_old+length(L)-1);
    temp_kulma = [temp_kulma; max(temp)-min(temp)];
end

Nappar_eff = Nappar;
date_tol_strict = 3;
if ( exist('dates_back','var') && length(dates_back) > 0 )
    dates_back = dates_back - dates_back(1);
    ang_back = rad2deg(ang_back);
end

% Phase correction
leader_phasecorr

% Compute eta for each apparition
for i=1:Nappar_eff
    % The start and end indices of the data
    ind = sum(pointsperapp(1:i-1))+1;
    inde = ind+pointsperapp(i)-1;
    
    if ( Nappar > Nappar_eff )
        % The case of combined apparitions
        L = L_back;
    else
        L = L_back(ind:inde);
    end

    % If we have over 'wanted' data points, compute eta and A
    if ( length(L) >= wanted )
        Npisteita = [Npisteita; length(L)];
        L2 = L.^2;
        eta = std(L2)/mean(L2);
        A = [A; sqrt( 1 - (1/(sqrt(8)*eta)+1/2)^(-1) )];
    end
    
end

% Remove complex and non-finite amplitudes
if ( ~isreal(A) || any(~isfinite(A)) )
    A( A ~= real(A) ) = [];
    A(~isfinite(A)) = [];
end