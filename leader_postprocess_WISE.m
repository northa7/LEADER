% Visual postprocessing ("deconvolution") for the
% reconstructed joint distribution f(p, beta)

disp('Smoothing the solution...')

% Do we allow spreading in the p-axis
% (could be ok for synthetic data,
% not recommended for real data)
if ( exist('allow_p_spread', 'var') && allow_p_spread )
    dampen = 0.1;
else
    dampen = 1;
end

[pind,bind]=find(W==max(max(W)));
W_after = W;
% Dampen the bins when moving away from the peak
for i=1:size(W,1)
    for j=1:size(W,2)
        W_after(i,j) = W(i,j)/...
            ( ( dampen*abs(pind-i)+abs(bind-j)+1 )^1 );
    end
end

% Move the values of p to the right by a constant (fixed) step
Pshift = 0.1;
PP = [P(1), min(P(2:end)+Pshift, 1)];
% Make sure to keep the end of PP increasing
ind = find(PP == 1);
if ( length(ind) > 1 )
    temp = PP(ind(1)-1);
    for i=1:length(ind)-1
        % Linear progression between [temp, 1]
        PP(ind(i)) = temp + i/length(ind)*(1-temp);
    end
end
% The shift in beta direction is hard to predict, so no shift is made
BB = BETA;

% Smoothed contour
figure
contourf(PP,BB,W_after')
colorbar
xlabel('p'), ylabel('\beta')
drawnow