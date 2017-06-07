% Read Pan-STARRS1 data files

% Load the PS1 database
if ( ~exist('lcg_files', 'var') )
    load lcg_files_PS1
end

% If the object we want to read hasn't been declared separately,
% choose a random object from the database
if ( ~exist('mones', 'var') || isempty(mones) )
    mones = ceil(length(lcg_files)*rand);
end
tiednimi = lcg_files{mones};
fid = fopen(tiednimi);
clear mones

% The first row has the number of lightcurves (not in PS1)
%temp = textscan(fid, '%f', 1, 'Delimiter', '\n');
%Nasteroids = temp{1};
Nasteroids = 1;
data=cell(1, Nasteroids);
% Go through each lightcurve
for i=1:Nasteroids
    % First, a comment line
    temp=textscan(fid, '%s', 1, 'Delimiter', '\n');
    % The number of measurements is stated in the comment line
    temp=textscan(temp{1}{1}, '%f');
    norows=temp{1}(1);
    % Collect the measurement data
    data{i}=textscan(fid, '%f %f %f %f %f %f %f %f', norows);
    % Skip the last line change
    temp = textscan(fid, '%c', 1, 'Delimiter', '\n');
end
fclose(fid);

% To avoid bias, choose one of the lightcurves randomly
mones2 = ceil(Nasteroids*rand);
data = data{mones2};

% Assort the data
dates = data{1};
L_big = data{2};
e_sun = [data{3}, data{4}, data{5}];
e_earth = [data{6}, data{7}, data{8}];
% Normalize
e_sun = normr(e_sun);
e_earth = normr(e_earth);

% Remove the measurements where the phase angle > ang_tol
ang_tol = deg2rad(20);
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
date_tol = 3;
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
date_tol_appar = 60;
if ( exist('dates_back','var') && length(dates_back) > 0 )
    dates_back = dates_back - dates_back(1);
    ang_back = rad2deg(ang_back);
end

% Phase correction
leader_phasecorr

% Compute eta for each set
for i=1:Nappar_eff
    % The start and end indices of the data
    ind = sum(pointsperapp(1:i-1))+1;
    inde = ind+pointsperapp(i)-1;
    
    if ( Nappar > Nappar_eff )
        % The case of combined sets
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