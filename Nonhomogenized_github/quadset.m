function [ w,p ] = quadset( dim,flag )
%Weights and points of std int,tri,pr, rect for Gauss quad in given dim
%   An entry in flag is only used to tell code nargin = 2 and thus we have
%   rectangles to deal with. Maple revealed that the K_bi matrix is degree
%   2 in both xi_1 and xi_2 respectively.  Thus, if use a 3 point
%   quadrature in those variables, then \int^1_0\int^(1-x)_0a(x)b(y)dydx
%   leads to a a(x)*deg 3 in x to be integrated, gives a final deg 5. 

% dim == 1: 
% This from Schumaker
%p1 =   1-[.966234757, .830604693, .619309593, ...
%         .380690407, .169395307, .033765243]';
%w1  = .5*[.171324492, .360761573, .467913935,...
%         .467913935, .360761573, .171324492]';     
% This from Taylor and COV from [0,1] back to [-1,1]
p1_3 = .5+.5*[-sqrt(.6);0;sqrt(.6)];
w1_3 = .5*[5/9;8/9;5/9];

p1_6 = .5+.5*[-.932469514203152;-.661209386466265;-.238619186083197;...
            .238619186083197;.661209386466265;.932469514203152];
w1_6 = .5*[.171324492379170;.360761573048139;.467913934572691;...
         .467913934572691;.360761573048139;.171324492379170];

% dim == 2: Fubini the 1d.  Leads to
%(1-xi_j)w_i*w_j*f(xi_i,(1-xi_i)*xi_j))
% First write as i,j arrays, ie w_j means constant with change in i
w2 = diag(1-p1_3)*(repmat(w1_3',3,1).*repmat(w1_3,1,3));

p2x = repmat(p1_3,1,3);
p2y = diag(1-p1_3)*repmat(p1_3',3,1);

% dim == 3: Fubini to 3d. Leads to
% (1-xi_i)w_i*w_j*w_k*f(xi_i,(1-xi_i)*xi_j,xi_k)
% Remember we use a 6 point quadrature in z. 
% First write as 3 x 3 x 6 array. Recall repmat
% drops the given array at each point of the
% block indices.
w3 = repmat(reshape((1-p1_3).*w1_3,3,1,1),1,3,6).*... i index
     repmat(reshape(w1_3,1,3,1),3,1,6).*... j index
     repmat(reshape(w1_6,1,1,6),3,3,1); %k index

p3x= repmat(reshape(p1_3,3,1,1),1,3,6);
p3y= repmat(reshape(1-p1_3,3,1,1),1,3,6).*...
     repmat(reshape(p1_3,1,3,1),3,1,6);
p3z= repmat(reshape(p1_6,1,1,6),3,3,1);

% flag:  Rectangle 2d quadrature
wr = repmat(w1_6',6,1).*repmat(w1_6,1,6);
prx= repmat(p1_6,1,6);
pry= repmat(p1_6',6,1);


if nargin == 1
    %It's ok to use : because all arrays match
    %according to whether the underlying (i,j,k)
    %index matches.
    if dim == 1
        w = w1_6;
        p = p1_6;
    elseif dim == 2
        w = w2(:);
        p = [p2x(:) p2y(:)];
    elseif dim == 3
        w = w3(:);
        p = [p3x(:),p3y(:),p3z(:)];
    else 
        error('Dim must be 1,2, or 3');
    end
else
      w = wr(:);
      p = [prx(:) pry(:)];
end

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

