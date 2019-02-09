function [ V,p_sample ] = ...
         check_discsample( p,tri,E_st,rinc,thetainc,scaling )
%Compute value of S^0_1 spline E_st at sample points
%   Detailed explanation goes here


%Initialize values needed for loop
p = p*scaling;
p_sample = check_genpts(rinc,thetainc)';
p_sample = p_sample*scaling;
n_sample = size(p_sample,2);
V = zeros(n_sample,1);
[~, ~, ~, pt_min, pt_max] = ... So to not redundantly loop generate lists
check_findtri(p,tri,p_sample(:,1),0);

%Evaluate function at all sample points.
for i=1:n_sample
    
    %Identify the triangles in which p_sample's belong
    [tfound, ~, ...
        ~] = check_findtri(p,tri,p_sample(:,i),0,pt_min,pt_max);
    nodes = tri(:,tfound); %Vertices of the discovered triangle
    
    %Extract the x-y coordinates of tri's nodes and 
    %get E_st's nodal values there
    P = [p(:,nodes(1)) p(:,nodes(2)) p(:,nodes(3))];
    nodals = [E_st(nodes(1)) E_st(nodes(2)) E_st(nodes(3))];
    
    %For the shape functions are varphi1 = 1-\xi1-\xi2, 
    %varphi2 = \xi1, and \varphi3 = \xi2. The equations are
    %\vec{x} = P*varphi 
    %spline(\xi) = nodals*varphi and
    %varphi = [1;0;0] +[-1 -1;1 0;0 1]*[\xi1,\xi2]
    %Thus, 
    
    P_LI = (P*[-1 -1;1 0;0 1])^(-1);
    xi   =  P_LI*(p_sample(:,i)-P*[1;0;0]);
    V(i) = nodals*(...
           [1;0;0]+[-1 -1;1 0;0 1]*xi...
           );
    
end
    
end

