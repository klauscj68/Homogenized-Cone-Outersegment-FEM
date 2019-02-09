function [ prfound,bbest,b,pt_min,pt_max ] =...
           check_3dfindpr( p,pr,tri,p_cross,...
                           pfind,tol,pt_min,pt_max)
%Find first prism in mesh that contains pfind. Dom is horizontal rad R_t.
%   TO save time, pr should only have the 6 rows for vertex indices.
%   tri=pr(1:3,:) should be prebuilt. p_cross = p(1:2,:) should be
%   prebuilt.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the data
R_b = 3;
R_t = 1;
H   = 15;
nu  = 1;
n_chambers = 100;
epsilon_0 = H/((1+nu)*(n_chambers-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_p  = size(p,2);
n_pr = size(pr,2);
%pr   = pr(1:6,:);

%Generate pt_min and pt_max if they are not given
if nargin == 6
    %pr is vertex indices which is why it's called
    xcoords= [p(1,pr(1,:));p(1,pr(2,:));p(1,pr(3,:));... p(1, means x-coord
              p(1,pr(4,:));p(1,pr(5,:));p(1,pr(6,:))];
    ycoords= [p(2,pr(1,:));p(2,pr(2,:));p(2,pr(3,:));... p(2, means y-coord
              p(2,pr(4,:));p(2,pr(5,:));p(2,pr(6,:))];
    zcoords= [p(3,pr(1,:));p(3,pr(2,:));p(3,pr(3,:));... p(3, means z-coord
              p(3,pr(4,:));p(3,pr(5,:));p(3,pr(6,:))];
    
    pt_min = [min(xcoords);min(ycoords);min(zcoords)];
    pt_max = [max(xcoords);max(ycoords);max(zcoords)];
    
end

%Find which prisms have their boxes containing the given point
xbds = (pfind(1)>=pt_min(1,:)-tol).*(pfind(1)<=pt_max(1,:)+tol);
ybds = (pfind(2)>=pt_min(2,:)-tol).*(pfind(2)<=pt_max(2,:)+tol);
zbds = (pfind(3)>=pt_min(3,:)-tol).*(pfind(3)<=pt_max(3,:)+tol);

bds = xbds.*(ybds.*zbds);
cand= find(bds);

if isempty(cand)
    error('Point not found in bounding rectangles. Adjust tolerance')
end

%Find first prism/best prism that does actually contain the point
%Because these prisms are tensor product, if prism contains the x,y
%value, we already checked that it contains the z.  We'd be done.
%disp('Make sure prism list is same triangle tensor stacked on itself')
%p_cross = p(1:2,:);
%tri = pr(1:3,:);

%Now find out which triangle we actually belong to.  Do this
%computation by checking the signs of its barycentric coordinates
%relative each one of candidate triangles.

ncand = size(cand,2);
i = 1; %Next two are to setup running while loop
bbest = -inf;

while (bbest < 0)&&(i<=ncand)
    b = [1 1 1;... %System for barycentric coordinates
         p_cross(1,tri(1:3,cand(i))');
         p_cross(2,tri(1:3,cand(i))')]^(-1)...
        *[1;pfind(1:2,1)];
    btest = b.*(b<=0);
    btest = sum(btest);
    if btest == 0 %All bary coordinates were positive.
        prfound = cand(i);
        bbest = 0;
    elseif bbest < btest %We have found a better approximate
                         %containing triangle
        prfound = cand(i);
        bbest = btest;
    end
    i = i+1;
end

end

