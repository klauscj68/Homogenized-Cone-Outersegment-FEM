function [ M ] = rect_M( X,Y,Z,weights,pts )
%Assemble mass matrix over surface rectangle element by Gauss Quadrature
%   X and Y are 4 by 1 column vectors that encode in first two entries
%   the coordinates of bottom segment's nodes and in final two
%   the top segment's nodes.  The nodes should be listed in either a
%   counter or counter-clockwise tour and not bottom two then top two with
%   vertical alignment.  That could lead to an isoparameteric mapping with
%   a cross.  The Z vector is a 2 by 1 column
%   vector which encodes the z-height of the bottom and top faces
%   respectively. Recall the reference rectangle is at 0, e_1, e_1+e_2 and
%   e_2.

% The maple file's output computes the 
% the \phi_{i}*\phi{j}
% products as a function of prism nodes
% and the standard prisms \xi coordinates.
% It also outputs the transformation's
% Jacobian squared up to possibly needing to
% take an absolute value sign. 

% 1d Weights for 6 point Gaussian quadrature
% and samples and is exact up to polys of deg
% 11
%w  = [.171324492, .360761573, .467913935,...
%      .467913935, .360761573, .171324492]';
%  
%xi = 1 - [.966234757, .830604693, .619309593, ...
%         .380690407, .169395307, .033765243]';
%
% Write \int_{I^2}f(x,y)dxdy = \int^1_0[.5*w_j*f(x,xi_{j})]dy
% = .25*w_i*w_j*f(xi_i,xi_j)
%
% Tensor product quadrature points and weights. The points are
% arranged in Matlab's indexing scheme so that the output feval
% vector matches going down the columns of weights.
% This makes it play with reshape.
%weights = .25*...
% [w(1)*w(1) w(1)*w(2) w(1)*w(3) w(1)*w(4) w(1)*w(5) w(1)*w(6);...
%  w(2)*w(1) w(2)*w(2) w(2)*w(3) w(2)*w(4) w(2)*w(5) w(2)*w(6);...
%  w(3)*w(1) w(3)*w(2) w(3)*w(3) w(3)*w(4) w(3)*w(5) w(3)*w(6);...
%  w(4)*w(1) w(4)*w(2) w(4)*w(3) w(4)*w(4) w(4)*w(5) w(4)*w(6);...    
%  w(5)*w(1) w(5)*w(2) w(5)*w(3) w(5)*w(4) w(5)*w(5) w(5)*w(6);...
%  w(6)*w(1) w(6)*w(2) w(6)*w(3) w(6)*w(4) w(6)*w(5) w(6)*w(6)];
%  
%pts = [...
%[xi(1);xi(1)] [xi(2);xi(1)] [xi(3);xi(1)] [xi(4);xi(1)] [xi(5);xi(1)] ...
%                                                        [xi(6);xi(1)]...
%[xi(1);xi(2)] [xi(2);xi(2)] [xi(3);xi(2)] [xi(4);xi(2)] [xi(5);xi(2)] ...
%                                                        [xi(6);xi(2)]...
%[xi(1);xi(3)] [xi(2);xi(3)] [xi(3);xi(3)] [xi(4);xi(3)] [xi(5);xi(3)] ...
%                                                        [xi(6);xi(3)]...
%[xi(1);xi(4)] [xi(2);xi(4)] [xi(3);xi(4)] [xi(4);xi(4)] [xi(5);xi(4)] ...
%                                                        [xi(6);xi(4)]...
%[xi(1);xi(5)] [xi(2);xi(5)] [xi(3);xi(5)] [xi(4);xi(5)] [xi(5);xi(5)] ...
%                                                        [xi(6);xi(5)]...
%[xi(1);xi(6)] [xi(2);xi(6)] [xi(3);xi(6)] [xi(4);xi(6)] [xi(5);xi(6)] ...
%                                                        [xi(6);xi(6)] ];
%
% The integrand to be taken over a standard 
% reference rectangle.
M = zeros(4,4);
I=eye(4);

% Turns out the Maple file was built so that 
% shape function coordinates matched rect verts
% [0,0], [0,1], [1,0], [1,1]
% whereas I've been assuming that they were
% listed as (counter-)clockwise, either way.
% I fix that now.

Perm = [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];
X = Perm*X;
Y = Perm*Y;


for i=1:4
    for j=1:4

% _preM doesn't need XX,YY,Z because on std element
% these guys pull back to interpolant functions on
% std element.  The geometry of physical elem is 
% encoded by the determinant.
f = @(xi1,xi2)  (I(i,:)*...
            (rect_nonhom_preM_bi([xi1;xi2])'*...
             rect_nonhom_preM_bi([xi1;xi2])*...
            I(:,j)))*...
           sqrt(abs(...                 %SA Jacobian
           rect_nonhom_predetJ_bi_patch(X,Y,Z,[xi1;xi2])...
           ));

eval = arrayfun(f,pts(:,1),pts(:,2));
M(i,j) = sum(weights.*eval);

    end
end

end

