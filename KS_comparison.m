% -------------------------------------------------------------------------
%
% Compare two distributions
%
% Assumptions: marginal DFs of p and beta have already
% been reconstructed for both populations.
%
% Inputs: two .mat files, which contain the following variable names:
% - in file01: P1, Pmargin1, BETA1, Bmargin1
% - in file02: P2, Pmargin2, BETA2, Bmargin2 
%
% Outputs:
% - DvalueP: the statistical differences for the CDFs of p (L1, L2, L_inf)
% - DvalueB: the statistical differences for the CDFs of beta (as above)
% - results: DFs and CDFs so that they can be plotted
%
% -------------------------------------------------------------------------


function [DvalueP, DvalueB, results] = KS_comparison(file01, file02)

load(file01)
load(file02)


% Normalize the marginal DFs
Pmargin1 = Pmargin1/sum(Pmargin1);
Pmargin2 = Pmargin2/sum(Pmargin2);
Bmargin1 = Bmargin1/sum(Bmargin1);
Bmargin2 = Bmargin2/sum(Bmargin2);


% Shape elongation p
% ------------------

% Add 0 and 1 in the ends to avoid NaNs
P1 = [0, P1, 1];
P2 = [0, P2, 1];
% CDFs (note: Pmargin is a column vector)
CP1 = [0, cumsum(Pmargin1')/sum(Pmargin1), 1];
CP2 = [0, cumsum(Pmargin2')/sum(Pmargin2), 1];
% Interpolate the CDF of P2 at the grid points of P1
CP2i = interp1(P2,CP2,P1);

% The statistical differences
DvalueP = [norm(CP1-CP2i, 1)/4, norm(CP1-CP2i), 2*norm(CP1-CP2i, inf)];


% Splin latitude beta
% -------------------

% Add 0 and pi/2 in the ends to avoid NaNs
BETA1 = [0, BETA1, pi/2];
BETA2 = [0, BETA2, pi/2];
% CDFs (note: Bmargin is a row vector)
CB1 = [0, cumsum(Bmargin1)/sum(Bmargin1), 1];
CB2 = [0, cumsum(Bmargin2)/sum(Bmargin2), 1];
% Interpolate the CDF of BETA2 at the grid points of BETA1
CB2i = interp1(BETA2,CB2,BETA1);

% The statistical differences
DvalueB = [norm(CB1-CB2i, 1)/4, norm(CB1-CB2i), 2*norm(CB1-CB2i, inf)];

% Plots
figure
subplot(2,1,1)
plot(P1, [0;Pmargin1;0], 'b-', 'LineWidth', 3)
hold on
plot(P2, [0;Pmargin2;0], 'r-.', 'LineWidth', 3)
xlabel('p')
ylabel('DF of p')
title(['Differences: D(L^1) = ', num2str(DvalueP(1)), ...
    ', D(L^2) = ', num2str(DvalueP(2)), ...
    ', D(L^{\infty}) = ', num2str(DvalueP(3))])
subplot(2,1,2)
plot(P1, CP1, 'b-', 'LineWidth', 3)
hold on
plot(P1, CP2i, 'r-.', 'LineWidth', 3)
xlabel('p')
ylabel('CDF of p')

figure
subplot(2,1,1)
plot(BETA1, [0,Bmargin1,0], 'b-', 'LineWidth', 3)
hold on
plot(BETA2, [0,Bmargin2,0], 'r-.', 'LineWidth', 3)
xlabel('\beta')
ylabel('DF of \beta')
title(['Differences: D(L^1) = ', num2str(DvalueB(1)), ...
    ', D(L^2) = ', num2str(DvalueB(2)), ...
    ', D(L^{\infty}) = ', num2str(DvalueB(3))])
subplot(2,1,2)
plot(BETA1, CB1, 'b-', 'LineWidth', 3)
hold on
plot(BETA1, CB2i, 'r-.', 'LineWidth', 3)
xlabel('\beta')
ylabel('CDF of \beta')

if ( nargout == 0 )
    disp('Statistical differences (format: [L1-norm, L2-norm, inf-norm]):')
    DvalueP
    DvalueB
    clear DvalueP DvalueB
end

% Collect the DFs and CDFs in a cell table
if ( nargout == 3 )
    results = cell(2,8);
    
    results{1,1} = {'p1', 'p2', 'DF of p1', 'DF of p2', ...
        'CDF of p1', 'CDF of p2', 'CDF of p2 (interp)'};
    results{1,2} = P1;
    results{1,3} = P2;
    results{1,4} = [0,Pmargin1',0];
    results{1,5} = [0,Pmargin2',0];
    results{1,6} = CP1;
    results{1,7} = CP2;
    results{1,8} = CP2i;
    
    results{2,1} = {'beta1', 'beta2', 'DF of beta1', 'DF of beta2', ...
        'CDF of beta1', 'CDF of beta2', 'CDF of beta2 (interp)'};
    results{2,2} = BETA1;
    results{2,3} = BETA2;
    results{2,4} = [0,Bmargin1,0];
    results{2,5} = [0,Bmargin2,0];
    results{2,6} = CB1;
    results{2,7} = CB2;
    results{2,8} = CB2i;
end