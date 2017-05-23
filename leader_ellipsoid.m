% Go through the faces of an asteroid model to
% compute shape elongation p = b/a

% Get the normals, areas and centers of triangulated faces
normaali = zeros(size(F));
ala = zeros(length(normaali), 1);

% Go through the faces
for i=1:size(F,1)
    p1 = R(F(i,1),:);
    p2 = R(F(i,2),:);
    p3 = R(F(i,3),:);
    
    a1 = p2-p1;
    a2 = p3-p2;
    
    % Cross product and determinant are actually faster to compute
    % without function calls
    %temp = cross2(a1,a2);
    temp = [a1(2).*a2(3)-a1(3).*a2(2), ...
        a1(3).*a2(1)-a1(1).*a2(3), ...
        a1(1).*a2(2)-a1(2).*a2(1)];
    normtemp = sqrt(temp(1)^2+temp(2)^2+temp(3)^2);
    normaali(i,:) = temp/normtemp;
    ala(i) = 1/2*normtemp;
end

% Grids for the angles
%epsilon = linspace(0, pi, 20);
%rotangles = linspace(0, 2*pi, 10);
%lambda = linspace(0, 2*pi, 15);
%mu = zeros(size(F,1),1);

% Get the semiaxes (based on the code written by Mikko Kaasalainen)
X = R(:,1);
Y = R(:,2);
Z = R(:,3);
phi = (1:180)/180*pi;

% x-direction
xphi = zeros(size(phi));
for i=1:length(phi)
    xx = X.*cos(phi(i)) + Y.*sin(phi(i));
    xphi(i) = max(xx)-min(xx);
end
a = max(xphi);
phimax = phi(a == xphi);

% y-direction
yy = Y.*cos(phimax)-X.*sin(phimax);
b = max(yy)-min(yy);

% z-direction
c = max(Z)-min(Z);

semiax = [a, b, c]/c;
%p = mean([b c])/a
p = b/a;