function [ F,p_sample ] = ...
         check_discfunc( p,tri,M,E_st,rinc,thetainc,scaling )
%Compute \int_{D}E_st(-) on dirac splines concentrated at samples' tri's  
%   p, tri are the mesh given as rows.  M is the mass matrix over the disc.  
%   E_st is the column vector whose coefficients are values of spline at
%   the corresponding index p-node. rinc, thetainc are the increments used
%   to generate the unit disc sample grid.  scaling is a factor that takes
%   you from a unit disc into the radius of the activation disc.

%When checking hom use the command line
%[F_hom,p_sample_hom] = check_discfunc(p_pd_sd(1:2,:)/R_b,t_pd_sd(1:3,:),M{1},E_st_sd{1}(i,:)',.05,pi/20,R_b*Z_scaling);
%Do this because the p_pd_sd(1:2,:) were meshed in disc radius R_b and need
%to bring it back to disc radius 1 so p_sample and nodes are scaled by same
%factor to get activation radius. The last factor is R_b*Z_scaling because
%it needs to a disc of radius 1 to a disc of radius R_b because Z_scaling
%takes a disc of radius R_b into the activation radius.

%When checking nonhom use the command line
%[F_nonhom,p_sample_nonhom] = check_discfunc(p_pd(1:2,:),t_pd(1:3,:),M{1},E_st_sd{1}(i,:)',.05,pi/20,z_scaling);
%Do this because p_pd and scaling points are already in unit disc and just
%need to be rescaled to activation site.


%Initialize values needed for loop
p = p*scaling;
p_sample = check_genpts(rinc,thetainc)';
p_sample = p_sample*scaling;
n_sample = size(p_sample,2);
F = zeros(n_sample,1);
[~, ~, ~, pt_min, pt_max] = ...
check_findtri(p,tri,p_sample(:,1),0);

for i=1:n_sample
    %Identify the triangles in which p_sample's belong
    [tfound, ~, ...
        ~] = check_findtri(p,tri,p_sample(:,i),0,pt_min,pt_max);
    
    %Compute L1 norm of spline whose only nonzero nodals 
    %are the vertices of triangle to which point belongs
    %at which spline is unity. 
    vert1 = tri(1,tfound);
    vert2 = tri(2,tfound);
    vert3 = tri(3,tfound);
    l_1 = sum(...
        [M(:,vert1) M(:,vert2) M(:,vert3)]...
        *[1;1;1]);
    
    %Integrate E_st against s/||s||_{1}. The actual
    %value is s^tr/||s||_{1}*M*E_st, which now use
    %that s is nonzero only at those 3 nodes and is
    %unity there mod normalization. (unity gives sum)
    ram = M*E_st;
    F(i) = ram(vert1) + ram(vert2) + ram(vert3);
    F(i) = F(i)/l_1;  
end

end

