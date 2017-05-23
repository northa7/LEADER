% Plots for the reconstructed DFs of p and beta

% Drawing the CDF of A and the function series sum(w_ij*F_ij)
% into the same plot to see how well the function series
% approximates the CDF
figure
plot(Asort, CDFA, 'bo')
hold on
plot(Asort, M*W_back, 'rx')
grid on
xlabel('A')
legend('CDF of A', 'sum of w_{ij} F_{ij}')
title(['Relative error of the fit: ', num2str(relerr)])
drawnow

% Plot the occupation numbers (weights) of the bins
figure
% W' is needed due to syntax
surf(P, BETA, W')
hold on
% Highlight the highest peak
plot3(pmax,betamax,max(max(W)),'rx', 'LineWidth',10)
% The actual location of the peak (if known)
if ( exist('p','var') && exist('beta','var') && ...
        length(p) == 1 && length(beta) == 1 )
    % W' is needed due to syntax
    w_interp = interp2(P_Gr, BETA_Gr, W', p, beta);
    plot3(p,beta, w_interp, 'mx', 'LineWidth',10)
end
set(gca, 'XLim', [0 1])
set(gca, 'YLim', [0 pi/2])
xlabel('p'), ylabel('\beta'), zlabel('w')
drawnow

% Contour plot
figure
contourf(P,BETA,W')
colorbar
xlabel('p'), ylabel('\beta')
drawnow

% Marginal DFs
if ( ~exist('marginplots','var') || isempty(marginplots) )
    marginplots = false;
end

Pmargin=sum(W,2);
Bmargin=sum(W,1);
if ( marginplots )
    weightbar(P, Pmargin, 0, 1, 10)
    xlabel('p')
    ylabel('w')
    % Scaling
    temp = max(get(gca, 'YLim'));
    set(gca, 'YLim', [0 max(temp, 0.2)])
    drawnow

    weightbar(BETA, Bmargin, 0, pi/2, 16)
    xlabel('\beta')
    ylabel('w')
    % Scaling
    temp = max(get(gca, 'YLim'));
    set(gca, 'YLim', [0 max(temp, 0.2)])
    drawnow
end