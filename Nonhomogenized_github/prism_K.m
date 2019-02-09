function [ K ] = prism_K( X,Y,Z,weights,pts )
%Assemble stiffness matrix over an element by Gauss quadrature
%   X and Y are 6 by 1 column vectors that encode in first three entries
%   the coordinates of bottom triangular face's nodes and in final three
%   the top face's nodes.  The i and i+3 positions should be bottom and top
%   corresponding nodes respectively. The Z vector is a 2 by 1 column
%   vector which encodes the z-height of the bottom and top faces
%   respectively. Recall the reference prism is cross section a right
%   triangle at e_{2}, 0, e_{1} and z-height 1.

% The maple file's output computes the 
% the \nabla_x\phi_{i}*\nabla_x\phi{j}
% products as a function of prism nodes
% and the standard prisms \xi coordinates.
% It also outputs the transformation's
% Jacobian.

% The integrand to be taken over a standard reference prism
K = zeros(6,6);
I = eye(6);

% By inspection of Maple, the isoparametric map has
% Jcb = Jcb(\xi_{3}) and of poly order 2. By the tensorial
% nature of basis functions the mass matrix bi*bj
% decomposes as \int_{RefPr}\varphi(\xi1,\xi2)\psi(\xi3)J(\xi3)dxi
% deg varphi = 2 deg psi = 2, deg J =2. We have
% \int_{RefPr}\varphi(\xi1,\xi2)\psi(\xi3)J(\xi3)dxi = 
% \int^1_0\psi J(xi3)\int_{Std}\varphi dxidxi3 = 
%[1/3*\Sum_{mid\in tri}\varphi(mid)]*
% 1/2*[5/9*\psiJ((1-sqrt(.6))/2) + 8/9\psiJ(1/2) +
%     5/9*\psiJ((1+sqrt(.6))/2)]
% The above is a COV to a [-1,1] 3 pt Gaussian quad rule

% Build a weights matrix to match the integrand evaluated at
% G_{1d} by G_{2d} point pairs.

%weights = [1/2*5/9*1/3 5/54 5/54;...
%           1/2*8/9*1/3 8/54 8/54;...
%           1/2*5/9*1/3 5/54 5/54];
%       
% These are the evaluation points of 
% ref prism but arranged like the 
% 'down the column' indexing of matlab
% wrt to weights. This way may use 
% reshape to match weights.
%pts = [ [.5;0;(1-sqrt(.6))/2] [.5;0;.5] [.5;0;(1+sqrt(.6))/2]...
%        [0;.5;(1-sqrt(.6))/2] [0;.5;.5] [0;.5;(1+sqrt(.6))/2]...
%        [.5;.5;(1-sqrt(.6))/2] [.5;.5;.5] [.5;.5;(1+sqrt(.6))/2]] ;

for i=1:6
    for j=1:6

%f = @(xi1,xi2,xi3) (I(i,:)*...
%           reshape(prism_nonhom_preK_bi_patch(X,Y,Z,[xi1;xi2;xi3]),3,6)'*...
%           reshape(prism_nonhom_preK_bi_patch(X,Y,Z,[xi1;xi2;xi3]),3,6)*...
%            I(:,j))*...
%            abs(prism_nonhom_predetJ_bi(X,Y,Z,[xi1;xi2;xi3]));

eval = zeros(54,1);
for r=1:54
eval(r) = (I(i,:)*...
           reshape(prism_nonhom_preK_bi_patch(X,Y,Z,[pts(r,1),pts(r,2),pts(r,3)]),3,6)'*...
           reshape(prism_nonhom_preK_bi_patch(X,Y,Z,[pts(r,1),pts(r,2),pts(r,3)]),3,6)*...
           I(:,j))*...
           abs(prism_nonhom_predetJ_bi(X,Y,Z,[pts(r,1),pts(r,2),pts(r,3)]));       
end
        
%eval    = arrayfun(f,pts(:,1),pts(:,2),pts(:,3));
K(i,j)  = sum(weights.*eval);

    end 
end

end

