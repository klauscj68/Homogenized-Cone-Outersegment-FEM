function [ tfound,bbest,b,pt_min,pt_max]...
          = check_findtri( p,tri,pfind,tol,pt_min,pt_max )
%Given nodes p and triangles t, give a triangle that contains pfind
%   Below code is vectorized in the style of the cone code. Tolerance
%   allows the bounding min-max rectangles to not exactly contain the
%   point. tfound is the index of the triangle to which we belong. 
%   If bbest == -||b_{-}||_{l1} then our triangle was a best approximation.
%   If bbest is 0, then we have truly found a triangle. vector b is the
%   barycentric coordinates of pfind in that triangle ASSUMING tfound
%   was genuine and not a best approximation. 

%Build pt_min and pt_max if not pregenerated
if nargin == 4 

%First to each triangle find the cartesian rectangle
%containing it.  First row is x-bds and 2nd row is y-bds.
%Recall for linear functions like coordinates, as 
%gradient doesn't vanish/is constant, function
%guaranteed to demonstrate its max and min at
%highest order boudary of boundary, eg triangles
%at their vertices

coord_x = [p(1,tri(1,:));... %x-coordinates of vertices of triangles
           p(1,tri(2,:));...
           p(1,tri(3,:))];
       
coord_y = [p(2,tri(1,:));... %x-coordinates of vertices of triangles
           p(2,tri(2,:));...
           p(2,tri(3,:))];       
       
pt_min = [...              %Minimum x-y coordinates for each triangle
          min(coord_x);...
          min(coord_y)...
         ];

pt_max = [...
          max(coord_x);...%Maximum x-y coordinates for each triangle
          max(coord_y)...
         ];

end
%Now find which rectangles each point belongs indexed by the
%triangle it boxes. The idea is that
%pfind\in tri\subseteq\Rightarrow pfind\in rectangle

cand_min = (pfind(1)>=pt_min(1,:)-tol).*(pfind(2)>=pt_min(2,:)-tol);
cand_max = (pfind(1)<=pt_max(1,:)+tol).*(pfind(2)<=pt_max(2,:)+tol);

cand = find(cand_min.*cand_max); %Find tri's whose rect's contain point

if isempty(cand)
    error('Point not found in bounding rectangles. Adjust tolerance')
end

%Now find out which triangle we actually belong to.  Do this
%computation by checking the signs of its barycentric coordinates
%relative each one of candidate triangles.

ncand = size(cand,2);
i = 1; %Next two are to setup running while loop
bbest = -inf;

while (bbest < 0)&&(i<=ncand)
    b = [1 1 1;... %System for barycentric coordinates
         p(1,tri(1:3,cand(i))');
         p(2,tri(1:3,cand(i))')]^(-1)...
        *[1;pfind];
    btest = b.*(b<=0);
    btest = sum(btest);
    if btest == 0 %All bary coordinates were positive.
        tfound = cand(i);
        bbest = 0;
        return
    elseif bbest < btest %We have found a better approximate
                         %containing triangle
        tfound = cand(i);
        bbest = btest;
    end
    i = i+1;
end

end

