% Computing the brightness functions and eta values
% (WISE, synthetic version)

% With a 75% chance, choose beta close to B_PEAK, otherwise random
if ( rand <= 0.75 )
    beta = 0;
    while ( beta <= 0 || beta >= pi/2 )
        beta = B_PEAK + 0.05*randn;
    end    
else
    beta = pi/2*rand;
end

% Choose random lambda from [0, 2*pi]
lambda = 2*pi*rand;
% Choose random rotation period from 3..12 hours (units: days)
Trot = 9/24*rand+3/24;

% Get geometries from a WISE asteroid
lcg_read_synth_WISE
L_big2 = zeros(size(L_big));
% Hapke parameters
Hapke_param = [0.63, 0.04, 1.4, -0.4];
Hapke_rough = 20;
% Compute the brightness values for the synthetic asteroid
for i=1:length(L_big)
    % Inspection angles and coordinate transformation
    E_inert = e_earth(i,:)';
    E0_inert = e_sun(i,:)';
    ROTM = transform_mat(0, 2*pi/Trot, dates(i), 0, beta, lambda);
    omega = ROTM*E_inert;
    omega0 = ROTM*E0_inert;
    
    mu = normaali*omega;
    mu0 = normaali*omega0;
    % For a convex shape, the face is visible if mu(k) > 0
    visible = (mu > 0);
    visible0 = (mu0 > 0);
    
%     % Hapke brightness values for each face
%     Ltemp = zeros(length(mu),1);
%     for j=1:length(mu)
%         Ltemp(j) = hapke_bright(omega', omega0', mu(j), mu0(j), ...
%             Hapke_param, Hapke_rough);
%     end
%     % Remove instabilities
%     Ltemp(~isfinite(Ltemp)) = 0;
%     Ltemp = visible.*visible0.*ala.*Ltemp;
    
    % L-S & L brightness values for each face
    Ltemp = visible.*visible0.*ala.*(mu.*mu0./(mu+mu0)+0.1*mu.*mu0);
    
    L_big2(i) = sum(Ltemp);  
end
% Add noise to L
L_big2 = L_big2 + 0.01*mean(L_big2)*randn(size(L_big2));

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
while ( i < length(L_big2) )
    L = L_big2(i);
    for j=i+1:length(L_big2)
        % Accept the data if the time span <= date_tol, quit otherwise
        if ( dates(j) - dates(i) <= date_tol )
            L = [L; L_big2(j)];
            % Close the loop in case it's the last measurement of the LC
            if ( j == length(L_big2) )
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
        devi = std(L2)/mean(L2);
        A = [A; sqrt( 1 - (1/(sqrt(8)*devi)+1/2)^(-1) )];
    end
    
end

% Remove complex and non-finite amplitudes
if ( ~isreal(A) || any(~isfinite(A)) )
    A( A ~= real(A) ) = [];
    A(~isfinite(A)) = [];
end