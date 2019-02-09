function [ Area ] = faces_area( rsmall,rbig,H,nu,eps0 )
%Compute the surface area of the discs
%   

%tangent of cone angle and also the lambda(z) scaling factor relative the
%small radius
tanalpha = (rbig-rsmall)/H;
lambda = @(x) 1+tanalpha/rsmall*x;

%Total number of discs, via n(nu*eps0+eps0) = H
n = floor(H/( (1+nu)*eps0 ));

%The Height positions denoting switch between interdiscal and disc
zeta = zeros(1,2*n);
for i=1:n
    zeta(1,2*i) = .5*nu*eps0 + i*eps0 + (i-1)*nu*eps0; %Top of a disc
    zeta(1,2*i-1) = .5*nu*eps0 + (i-1)*eps0 + (i-1)*nu*eps0; %Bottom of a disc
end

%Compute the cone radii at specified heights
scaling = arrayfun(lambda,zeta);
Radius_Height = scaling.*(rsmall+zeros(1,2*n));

%Add up the area of all the faces
Area = 0;
for i=1:2*n
    Area = Area + pi*Radius_Height(1,i)^(2);
end

end

