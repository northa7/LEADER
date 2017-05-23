function [xmean, sumw] = weightbar(x, w, xmin, xmax, Np)

% Draws statistical bars of points x with weights w, with Np bars
% Outputs: mean of each bar, sum weight of each bar

if ( nargin < 4 )
    xmin = 0;
    xmax = 1;
    Np = 10;
elseif ( nargin < 5 )
    Np = 10;
end

xmean = zeros(Np,1);
sumw = zeros(Np,1);
% the length of each bar, limits for each bar
barlength = (xmax-xmin)/Np;
barmin = xmin+(0:Np-1)*barlength;
barmax = xmin+(1:Np)*barlength;

for i=1:Np
    xmean(i) = 1/2*(barmin(i)+barmax(i));
    if ( i < Np )
        barind = find( x >= barmin(i) & x < barmax(i) );
        sumw(i) = sum(w(barind));
    else
        barind = find( x >= barmin(i) & x <= barmax(i) );
        sumw(i) = sum(w(barind));
    end
end
figure
bar(xmean, sumw)

if ( nargout == 0 )
    clear xmean sumw
end