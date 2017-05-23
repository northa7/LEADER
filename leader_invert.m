% Statistical inversion for the parameters p and beta

NP = 20;
NBETA = 29;
% Generating bins: divide the set [0, 1] into NP equally spaced bins,
% and [0, pi/2] into NBETA equally spaced bins, and randomly generate
% a representative value for each bin
P = linspace(0.025, 0.975, NP);
BETA = linspace(0.025, 1.55, NBETA);

% Perturb P and BETA with a "truncated Gaussian distribution"
% to randomize the bins
%
% Experiment: p from bins of different sizes
% ---------------------------------------------
if ( exist('gridtype', 'var') && isequal(gridtype, 'dynamic') )
    temp = zeros(1,2);
    temp(1) = P(randi([1 5]));      % p < 0.25
    temp(2) = P(randi([6 8]));      % 0.25 < p < 0.4
    temp = [temp, P(9:end)];        % other ratios keep their scaling
    P = temp;
end
% ---------------------------------------------
coeff = 0.015;
temp = coeff*randn(size(P));
for i=1:length(temp)
    while ( abs(temp(i)) > 0.025 )
        temp(i) = coeff*randn;
    end
end
P = P + temp;
% Experiment: beta from bins of different sizes
% ---------------------------------------------
if ( exist('gridtype', 'var') && isequal(gridtype, 'dynamic') )
    temp = zeros(1,5);
    temp(1) = BETA(randi([1 6]));       % angles 0-18.5 (deg)
    temp(2) = BETA(randi([7 8]));       % angles 18.5-24.7
    temp(3) = BETA(randi([9 10]));      % angles 24.7-31
    temp(4) = BETA(randi([11 12]));     % angles 31-37.2
    temp(5) = BETA(randi([13 14]));     % angles 37.2-43.4
    temp = [temp, BETA(15:end)];        % other angles keep their scaling
    BETA = temp;
end
% ---------------------------------------------
temp = coeff*randn(size(BETA));
for i=1:length(temp)
    while ( abs(temp(i)) > 0.025 )
        temp(i) = coeff*randn;
    end
end
% Make sure the last BETA value doesn't exceed pi/2
BETA = min(BETA + temp, pi/2-eps);

% The data matrix M for the linear system M*w = C
M = zeros(length(Asort), length(P)*length(BETA));
ind = 1;
for j=1:length(P)
    for k=1:length(BETA)
        for i=1:length(Asort)
            if ( Asort(i) <= P(j) )
                M(i,ind) = 0;
            elseif ( Asort(i) < ...
                    sqrt( sin(BETA(k))^2+P(j)^2*cos(BETA(k))^2 ) )                
                M(i,ind) = pi/2 - acos( ...
                    sqrt( (Asort(i)^2-P(j)^2)/(1-P(j)^2) ) / ...
                    sin(BETA(k)) );
            else
                M(i,ind) = pi/2;
            end
        end
        ind = ind + 1;
    end
end

% Matrices RP and RB for regularizing
NN = length(P);
MM = length(BETA);
RP = zeros((NN-1)*MM, NN*MM);
RB = zeros(NN*(MM-1), NN*MM);

for k=1:size(RP,1)
    pindex = ceil(k/MM);
    RP(k,k)   = -1/(P(pindex+1)-P(pindex));        
    RP(k,k+MM) = 1/(P(pindex+1)-P(pindex));
end

for k=1:size(RB,1)
    ll = k + ceil(k/(MM-1)) - 1;
    bindex = mod(k, MM-1);
    if ( bindex == 0 )
        bindex = MM-1;
    end
    RB(k,ll)  = -1/(BETA(bindex+1)-BETA(bindex));
    RB(k,ll+1) = 1/(BETA(bindex+1)-BETA(bindex));
end

if ( ~exist('deltaP', 'var') || isempty(deltaP) )
    deltaP = 0.1;
end
if ( ~exist('deltaB', 'var') || isempty(deltaB) )
    deltaB = 1;
end
Mtilde = [M; sqrt(deltaP)*RP; sqrt(deltaB)*RB];
Ctilde = [CDFA; zeros(size(RP,1), 1); zeros(size(RB,1), 1)];

% Solve the extended linear system Mtilde*w = Ctilde
% using a least squares method with a positivity constraint
disp('Solving the weights w_ij for the bins (p_i, beta_j)...')
W = lsqnonneg(Mtilde, Ctilde);
disp('Solution obtained!')

% Reshape the solution vector into a matrix
W_back = W;
W = reshape(W, length(BETA), length(P))';
% If you wanted to reshape back into a vector:
% W = reshape(W', length(BETA)*length(P), 1);

% Find the peak
[pmax,betamax]=find(W==max(max(W)));
pmax = P(pmax);
betamax = BETA(betamax);

[P_Gr, BETA_Gr] = meshgrid(P, BETA);

relerr = norm(M*W_back-CDFA)/norm(CDFA);
disp(['The highest peak: P=', num2str(pmax), ', BETA=', num2str(betamax)])
disp(['Relative error: ', num2str(relerr)])