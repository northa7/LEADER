% -------------------------------------------------------------------------
% Read asteroid vertices and facets from a file.
% The file name is given as a string input 'tiedosto'.
% 'plotflag' is an optional input about whether to
% draw a 3D surface plot of the asteroid.
% -------------------------------------------------------------------------

function [x,y,z,F] = damit_model(tiedosto, plotflag)

if ( ~exist('plotflag','var') || isempty(plotflag) )
    plotflag = true;
end

% Read the object file
fid = fopen(tiedosto);
% Row syntax: [ v/f, x, y, z ]
data = textscan(fid, '%s %f %f %f');
fclose(fid);

% Modify the table so that it is easier to handle
data_vert = [data{2} data{3} data{4}];
ind1=1;
ind2=1;
% Find the number of vertices and faces for pre-allocation
for i=1:length(data{1})
    if ( data{1}{i} == 'f' )
        numv = i-1;
        numf = length(data{1})-numv;
        break
    end
end
V = zeros(numv,3);
F = zeros(numf,3);
for i=1:length(data{1})
    % 'v' in the front means it's a vertex
    if ( data{1}{i} == 'v')
        V(ind1,:)=data_vert(i,:);
        ind1=ind1+1;
    % 'f' in the front means it's a face
    elseif( data{1}{i} == 'f' )
        F(ind2,:)=data_vert(i,:);
        ind2=ind2+1;
    end
end
V=V/mean(mean(abs(V)));
x=V(:,1);
y=V(:,2);
z=V(:,3);

% Coordinate transformation from Cartesian to spherical
%r=sqrt(x.^2+y.^2+z.^2);
%t=atan2(y,x)';
%phi=acos(z./r)';

% 3D surface plot of the asteroid (optional)
if ( plotflag )
    figure
    trisurf(F, x, y, z)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    drawnow
end