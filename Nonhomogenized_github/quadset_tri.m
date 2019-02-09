function [ w,p ] = quadset_tri( n )
%For n_pts, get w's & pts for gaussian quad on triangle
%   The expression \sum wi*f(pi) has 3n unknowns, so expect it can
%   integrate polynomials up to degree d where (d+1)^2 <= 3*n. Note that
%   the degree d is x and y jointly.

d = floor(sqrt(3*n)-1);

% There are 3*n-(d+1)^2 remaining degrees of freedom,
% that need to be filled, do this in the slots about center
% of infty degree is leq d+1 is up to
% (d+2)^2-(d+1)^2 = 2*(d+1)+1

% Setup I and J
I = [];
for i=0:d
    I = [I; i+zeros(d+1,1)];
end

J = cumsum(ones(d+1,1))-1;
J = repmat(J,d+1,1);

% Now fill in remaining 3*n-(d+1)^2 slots
ram_I = [0:d+1 (d+1)*ones(1,d+1)]';
ram_J = [(d+1)*ones(1,d+1) d+1 (d+1)-(1:d+1)]';

rest = 3*n-(d+1)^2;
if rest ~= 0
    I = [I;ram_I( (1:rest)' )];
    J = [J;ram_J( (1:rest)' )];
end

if (size(I,1)~= 3*n) || (size(I,1)~=size(J,1))
    error('The basis elems are not dim of constraints.');
end

% These values quadrature should hit.
Bvals = zeros(3*n,1);
for k=1:3*n
    Bvals(k) = basis_int(I(k),J(k));
end

% Setup initial guesses for weights and pts
w0 = rand(n,1);
p = zeros(n,2);
for i=1:n
    p(i,1) = rand(1);
    p(i,2) = p(i,1)*rand(1);
end

X0  = [w0 p];
Fun = @(X) F(X(:,1),X(:,2:3),I,J,Bvals);

V = fsolve(Fun,X0);

w = V(:,1);
p = V(:,2:3);


end

function [int] = basis_int(i,j)
    % Integrate x^i*y^j*z^k over the standard prism 
    % Handled by integrating z and x-y vars separately

    prod_bx = 1;
    if (i ~= 0) && (j ~= 0)
        %\int_{Tri} x^i*y^j = \int^1_0\int^(1-x)_0 x^i*y^j dydx
        %= 1/(j+1)\int^1_0 x^i*(1-x)^(j+1)dx
        %= 1/(j+1)\int^1_0 \sum^(j+1)_(k=0)nCr(j+1,k)(-1)^k*x^(i+k) dx
        %= 1/(j+1)\sum^(j+1)_(k=0) (-1)^k/(k+i+1)*nCr(j+1,k)
        sum = 0;
        for sk = 0:j+1
            sum = sum + ((-1)^sk)/(sk+i+1)*nchoosek(j+1,sk);
        end
        prod_bx = 1/(j+1)*sum;
        
    elseif (i~=0) && (j==0)
        %\int_{Tri} (x^i)dx = \int^1_0*\int^(1-y)_0 x^i*dxdy
        %  = 1/(i+1)*\int^1_0 (1-y)^(i+1)dy = 1/(i+1)*1/(i+2)
        prod_bx = 1/( (i+1)*(i+2) );
        
    elseif (i==0) && (j~=0)
        prod_bx = 1/( (j+1)*(j+2) );
        
    elseif (i==0) && (j==0)
        %\int_{Tri} 1 dxdy = .5*base*height
        prod_bx = .5;
    end

    int = prod_bx;
end

function [V] = F(weights,pts,I,J,Bints)
    % Compute the value of the system to be set to 0
    % by satisfactory Gaussian weights and points
    % that are to achieve the integrals of the 
    % given [I,J,K] basis elem: x^i*y^j*z^k
    % [I,J,K] are nd x 1 col vectors.
    % weights is a n x 1 col vector as pts is a n x 3
    % vector of points on std prism. 
    
    nd = size(I,1);
    V = zeros(nd,1);
    for i=1:nd
        V(i) = sum(...
                      weights.*(...
                         pts(:,1).^(I(nd)).*...
                         pts(:,2).^(J(nd))...
                  ))...
               - Bints(i);
    end
end

