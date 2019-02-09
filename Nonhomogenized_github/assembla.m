% costruzione delle matrici di massa e rigidezza di ciascuna regione: vol, sd, os, in
function [M_vol, K_vol, M_fo,M_hd, M_lc, M_pan,M_sl, M_ch,...
          M_gl, K_gl, Sigma_cone,Sigma_sl,Sigma_hd, Volume_cone] =...
    assembla( p_3d,f_3d,pr_3d )
    %,n_pd, n_tri, n_sl, n_fo, n_sez, n_sd, ...
    %p_pd, t_pd, Z_s, sl2pd, fo2pd, sd2sez, z_scaling )
% DOUBLE CHECK .5*nu*epsilon_0 vol2surf conversion is
% done.

% Read data
[ R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,sigma,n_chambers,...
           dz_0,n_sd,Z_sd,flag_geom_sp, ...
           taglia,n_sez,tol_R,tol_angle,n_ref_cyto,n_ref_id,...
           metodo_cyto,...
           flag_model,flag_model_disc,...
           cc_R_st, D_R_st, cc_G_st, D_G_st, cc_E_st, D_E_st,...
           n_step_R,mu,lambda,nu_RE,nu_RG,...
           k_E, k_GE, PDE_s, mode_time, mode_space, n_sample,...
           n_Phi, Phi,...
           norma_inf,...
           peak_delta,...
           plot_niagara, plot_pool, inspect,...
           t_fin, n_step_t, metodo_id, solver, downsample,...
           u_tent,v_tent,tol_stat,...
           theta, alpha, tol_fix,...
           beta_dark,k_st, B_ca, F, f_ca,...
           cc_u,kk_u,cc_v,kk_v,...
           alpha_max,alpha_min,m_cyc,k_cyc,...
           j_cg_max,m_cg,K_cg,...
           j_ex_sat,K_ex] = data;

% Sigma_cone refers to the area
% over which the channels are distributed, eg HD's, LC's, Sliver, but NOT
% panels. 
flag_assemble = false;
if flag_assemble == true
% The equations for this assembly becomes first for cGMP:
% cc_u*(\int_{\Omega}\varphi_{i}\varphi_{j})u'_{j} + 
% k_u*(\int_{\Omega}\nabla\varphi_{i}\nabla\varphi_{j})u_{j}
% = .5*nu_RE*epsilon_0*(-beta_dark\int_{dhd}\varphi_{i}\varphi_{j})u_{j}
% +(-k*_{\sigma,hyd}\int_{dactivation}\varphi_{i}\varphi_{j})(E_st.*u)_{j}+ %This one is [PDE*}_{\sigma} on activated disk
% (.5*nu_RE*epsilon_0)*(\int_{dhd}\varphi_{i}\alpha(v))
% Giovanni says above should actually be nu not nu_RE

% And second for Ca
% cc_v*(\int_{\Omega}\varphi_{i}\varphi_{j})v'_{j} +
% k_v*(\int_{\Omega}\nabla\varphi_{i}\nabla\varphi_{j})v_{j}
% = (-1/(B_CA*F)\int_{sl}\varphi_{i}J_{ex}(v)) 
% + (1/(2*B_CA*F)\int_{sl}f_{CA}\varphi_{i}J_{cG}(u)

% Rewritten in terms of the mass and stiffness matrices this
% cc_u*M_vol[u'] + 
% k_u*K_vol[u] +.5*beta_dark*nu_RE*epsilon_0*M_fo[u] +...
% k*_{sigma,hyd}M_{hd}[E_st.*][u]=...
% .5*nu_RE*epsilon_0*M_fo[\alpha(v)]

% cc_v*M_vol[v'] + k_v*K_vol[v] =
% -1/(B_Ca*F)(M_fo+M_sl)[J_ex(v)] + f_{Ca}/(2*B_Ca*F)(M_fo+M_sl)[J_cG(u)]

% Together as one vector system, this reads
% [cc_u*M_vol zeros; zeros cc_v*M_vol][u';v'] +
% [k_u*K_vol + .5*beta_dark*nu_RE*epsilon_0*M_fo+k*_{sigma,hyd}M_{sd}[E_st.*] zeros;...
%  zeros k_v*K_vol][u;v] = 
%  [.5*nu_RE*epsilon_0*M_fo[\alpha(v)];...
%  -1/(B_Ca*F)(M_fo+M_sl)[J_ex(v)] + f_{Ca}/(2*B_Ca*F)(M_fo+M_sl)[J_cG(u)]]

% output

% M_vol, M_sd and M_os, are the mass matrices in the volume, sliver and special disc
% note: these matrices are not multiplied by any coefficient
% they are needed both to create the global matrices, 
% and the loading vector,
% the loading vector in each region is obtained as the product between the mass matrix relevant to the region and the nodal vector of the load 

fprintf('\nAssemble mass and stiffness matrices of vol, sd, sl\n')

% Number of volume degrees of freedom
n_p3d = size(p_3d,2);

% number of prisms
n_pr = size(pr_3d,2);

% volume
% initialization
M_vol=zeros(36*n_pr,3);   % There are 6 by 6 basis node                                              % Is this the right number?
K_vol=zeros(36*n_pr,3);   % pairs over std prism.
Volume_cone = 0;
% number of nonzero elements
nnz_K=0;
nnz_M=0;

% Quadrature data for reference element
% weights = [1/2*5/9*1/3 5/54 5/54;...
%          1/2*8/9*1/3 8/54 8/54;...
%           1/2*5/9*1/3 5/54 5/54];
       
%pts = [ [.5;0;(1-sqrt(.6))/2] [.5;0;.5] [.5;0;(1+sqrt(.6))/2]...
%        [0;.5;(1-sqrt(.6))/2] [0;.5;.5] [0;.5;(1+sqrt(.6))/2]...
%        [.5;.5;(1-sqrt(.6))/2] [.5;.5;.5] [.5;.5;(1+sqrt(.6))/2]] ;

[weights,pts] = quadset(3);
% prisms in the volume
    for t=1:n_pr
        % indices of the triangular node in the pd numbering and of the
        % segment in the z-direction mesh
        disp(['Progress through volume elements: ' num2str(t/n_pr*100)]); 
        ind_vol=pr_3d(1:6,t);                                                %Would automatically be ind_vol too.%
        % volume dof, in the progressive numbering of the volume             % local dof in the prism element?  
        % nodal coordinates in the cylindrical domain (will be rescaled to   % Point is to match the ref prism nodes 
        % get the coordinates in the cone domain)                            % to right index in physical space: ind_vol
        %X=p_pd(1,ind_tri); Not needed because now nodes rescaled            %Only question is is sort important 
        %Y=p_pd(2,ind_tri); By Matlab!
        %Z=Z_s(ind_sez);
        %lambda=z_scaling(ind_sez);
        X=p_3d(1,ind_vol);
        Y=p_3d(2,ind_vol);
        Z=p_3d(3,ind_vol([1 4]))';  %Grab the two distinct z heights
        % element matrices
        [K_elem]=prism_K(X',Y',Z',weights,pts);                                     % Replace with new routine
        [M_elem]=prism_M(X',Y',Z',weights,pts);
        % Aggregate cytosolic volume
        Volume_cone = Volume_cone + abs(ones(1,6)*M_elem*ones(6,1));
        %finds the nonzero elements of K
        [I,J,V]=find(K_elem);
        nuovi_nnz=length(I);
        % stacks in K_vol
        K_vol(nnz_K+1:nnz_K+nuovi_nnz,:)=[ind_vol(I), ind_vol(J), V];
        nnz_K=nnz_K+nuovi_nnz;
        % finds the nonzero elements of M
        [I,J,V]=find(M_elem);
        nuovi_nnz=length(I);
        % stacks in M_vol
        M_vol(nnz_M+1:nnz_M+nuovi_nnz,:)=[ind_vol(I), ind_vol(J), V];
        nnz_M=nnz_M+nuovi_nnz;
    end

% create M_vol and K_vol by adding contributes relevant to the same node     %They just mean you only have 3 lists until sparse makes matrix%
M_vol=sparse(M_vol(1:nnz_M,1),M_vol(1:nnz_M,2),M_vol(1:nnz_M,3),n_p3d,n_p3d);%And adds together entries belonging to same node pair%
K_vol=sparse(K_vol(1:nnz_K,1),K_vol(1:nnz_K,2),K_vol(1:nnz_K,3),n_p3d,n_p3d);
disp('M_vol and K_vol are assembled');


%%%%%%
% Identify horizontal faces by f(5,:)==-1 
% for bottom and top. Identify
% chamber patch exteriors by f(4,:)~=0 and 
% f(5,:)==1 union f(5,:)==-2. Identify panels by
% f(5,:) == 3||4. Identify sliver by f(5,:)==5.
%%%%%%

%Initialize Sigma_cone
Sigma_cone = 0;
Sigma_hd   = 0;                                                    
% Mass matrix for horizontal disc boundaries
% Find bd faces and initialize M
found_bds = find(f_3d(5,:)==-1);
n_bds     = size(found_bds,2);

M_hd = zeros(9*n_bds,3);  %These are triangles
nnz_M= 0;
    for t=1:n_bds
        disp(['Progress through horizontal boundaries: ' num2str(t/n_bds*100)]);
        % nodal indices of the triangle (in the ordering of the pd mesh)
        ind_tri=f_3d(1:3,found_bds(t));
        % x-y nodal coordinates in cone OUTPUT AS ROWS!
        X=p_3d(1,ind_tri);
        Y=p_3d(2,ind_tri);
        % element matrices
        %[K_elem]=tri_K(X,Y);
        [M_elem]=tri_M(X,Y);
        tri_area = [1 1 1]*M_elem*[1;1;1];
        tri_area = abs(tri_area);
        Sigma_cone = Sigma_cone+tri_area;
        Sigma_hd   = Sigma_hd + tri_area;
        % finds nonzero elements of K
        %[I,J,V]=find(K_elem);
        %nuovi_nnz=length(I);
        % stacks in K
        %K(nnz_K+1:nnz_K+nuovi_nnz,:)=[ind_tri(I), ind_tri(J), V];
        %nnz_K=nnz_K+nuovi_nnz;
        % finds nonzero elements of M
        [I,J,V]=find(M_elem);
        nuovi_nnz=length(I);
        % stacks in M
        M_hd(nnz_M+1:nnz_M+nuovi_nnz,:)=[ind_tri(I), ind_tri(J), V];
        nnz_M=nnz_M+nuovi_nnz;
    end
    % somma i chiummi
    M_hd=sparse(M_hd(1:nnz_M,1),M_hd(1:nnz_M,2),M_hd(1:nnz_M,3),n_p3d,n_p3d);
    %K_sd{d}=sparse(K(1:nnz_K,1),K(1:nnz_K,2),K(1:nnz_K,3),n_pd,n_pd);
    disp('M_hd is assembled');

% Quadrature data for rectangular ref elem.
% 1d Weights for 6 point Gaussian quadrature
% and samples and is exact up to polys of deg
% 11
%w  = [.171324492, .360761573, .467913935,...
%      .467913935, .360761573, .171324492]';
%  
%xi = 1 - [.966234757, .830604693, .619309593, ...
%         .380690407, .169395307, .033765243]';
%
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

[weights,pts] = quadset(2,true);

% Mass matrix for lateral cone chamber boundary
% Find bd faces and initialize M
ram       = (f_3d(4,:)~=0).*(f_3d(5,:)==1); %chamber to cytosol
ram       = ram + (f_3d(5,:)==-2);          %sliver to cytosol
found_bds = find(ram);
clear ram
n_bds     = size(found_bds,2);

M_lc = zeros(16*n_bds,3);  %These are rectangles
nnz_M= 0;
    for t=1:n_bds
        disp(['Progress through lateral boundaries: ' num2str(t/n_bds*100)]);
        % nodal indices of the triangle (in the ordering of the pd mesh)
        ind_rect=f_3d(1:4,found_bds(t));
        % x-y nodal coordinates in cone OUTPUT AS ROWS!
        X=p_3d(1,ind_rect);
        Y=p_3d(2,ind_rect);
        Z=p_3d(3,ind_rect([1 3])); % Grab distinct z-coords
                                 % recall toured as cycle
        % element matrices
        %[K_elem]=tri_K(X,Y);
        [M_elem]=rect_M(X',Y',Z',weights,pts);
        rect_area  = [1 1 1 1]*M_elem*[1;1;1;1];
        rect_area  = abs(rect_area);
        Sigma_cone = Sigma_cone+rect_area;
        % finds nonzero elements of K
        %[I,J,V]=find(K_elem);
        %nuovi_nnz=length(I);
        % stacks in K
        %K(nnz_K+1:nnz_K+nuovi_nnz,:)=[ind_tri(I), ind_tri(J), V];
        %nnz_K=nnz_K+nuovi_nnz;
        % finds nonzero elements of M
        [I,J,V]=find(M_elem);
        nuovi_nnz=length(I);
        % stacks in M
        M_lc(nnz_M+1:nnz_M+nuovi_nnz,:)=[ind_rect(I), ind_rect(J), V];
        nnz_M=nnz_M+nuovi_nnz;
    end
    % somma i chiummi
    M_lc=sparse(M_lc(1:nnz_M,1),M_lc(1:nnz_M,2),M_lc(1:nnz_M,3),n_p3d,n_p3d);
    %K_sd{d}=sparse(K(1:nnz_K,1),K(1:nnz_K,2),K(1:nnz_K,3),n_pd,n_pd);
    disp('M_lc is assembled');

% Mass matrix for panels
% Find bd faces and initialize M
ram       = (f_3d(5,:)==3)+(f_3d(5,:)==4);
found_bds = find(ram);
clear ram
n_bds     = size(found_bds,2);

M_pan = zeros(16*n_bds,3);  %These are rectangles
nnz_M= 0;
    for t=1:n_bds
        disp(['Progress through panel boundaries: ' num2str(t/n_bds*100)]);
        % nodal indices of the triangle (in the ordering of the pd mesh)
        ind_rect=f_3d(1:4,found_bds(t));
        % x-y nodal coordinates in cone OUTPUT AS ROWS!
        X=p_3d(1,ind_rect);
        Y=p_3d(2,ind_rect);
        Z=p_3d(3,ind_rect([1 3])); % Grab distinct z-coords
                                 % recall toured as cycle
        % element matrices
        %[K_elem]=tri_K(X,Y);
        [M_elem]=rect_M(X',Y',Z',weights,pts);
        %X=p_3d(1,ind_rect);
        %Y=p_3d(2,ind_rect);
        % element matrices
        %[K_elem]=tri_K(X,Y);
        %[M_elem]=rectangule_M_gen(X,Y);
        %Sigma_cone = Sigma_cone+rect_area(X,Y);
        % finds nonzero elements of K
        %[I,J,V]=find(K_elem);
        %nuovi_nnz=length(I);
        % stacks in K
        %K(nnz_K+1:nnz_K+nuovi_nnz,:)=[ind_tri(I), ind_tri(J), V];
        %nnz_K=nnz_K+nuovi_nnz;
        % finds nonzero elements of M
        [I,J,V]=find(M_elem);
        nuovi_nnz=length(I);
        % stacks in M
        M_pan(nnz_M+1:nnz_M+nuovi_nnz,:)=[ind_rect(I), ind_rect(J), V];
        nnz_M=nnz_M+nuovi_nnz;
    end
    % somma i chiummi
    M_pan=sparse(M_pan(1:nnz_M,1),M_pan(1:nnz_M,2),M_pan(1:nnz_M,3),n_p3d,n_p3d);
    %K_sd{d}=sparse(K(1:nnz_K,1),K(1:nnz_K,2),K(1:nnz_K,3),n_pd,n_pd);  
    disp('M_pan is assembled');
    
% Mass matrix for outersliver
% Find bd faces and initialize M
ram       = (f_3d(5,:)==5);
found_bds = find(ram);
clear ram
n_bds     = size(found_bds,2);

Sigma_sl = 0;
M_sl = zeros(16*n_bds,3);  %These are rectangles
nnz_M= 0;
    for t=1:n_bds
        disp(['Progress through sliver boundaries: ' num2str(t/n_bds*100)]);
        % nodal indices of the triangle (in the ordering of the pd mesh)
        ind_rect=f_3d(1:4,found_bds(t));
        % x-y nodal coordinates in cone OUTPUT AS ROWS!
        X=p_3d(1,ind_rect);
        Y=p_3d(2,ind_rect);
        Z=p_3d(3,ind_rect([1 3])); % Grab distinct z-coords
                                 % recall toured as cycle
        % element matrices
        [M_elem]=rect_M(X',Y',Z',weights,pts);
        rect_area  = [1 1 1 1]*M_elem*[1;1;1;1];
        rect_area  = abs(rect_area);
        Sigma_cone = Sigma_cone+rect_area;
        Sigma_sl  =Sigma_sl + rect_area;
        % finds nonzero elements of K
        %[I,J,V]=find(K_elem);
        %nuovi_nnz=length(I);
        % stacks in K
        %K(nnz_K+1:nnz_K+nuovi_nnz,:)=[ind_tri(I), ind_tri(J), V];
        %nnz_K=nnz_K+nuovi_nnz;
        % finds nonzero elements of M
        [I,J,V]=find(M_elem);
        nuovi_nnz=length(I);
        % stacks in M
        M_sl(nnz_M+1:nnz_M+nuovi_nnz,:)=[ind_rect(I), ind_rect(J), V];
        nnz_M=nnz_M+nuovi_nnz;
    end
    % somma i chiummi
    M_sl=sparse(M_sl(1:nnz_M,1),M_sl(1:nnz_M,2),M_sl(1:nnz_M,3),n_p3d,n_p3d);
    %K_sd{d}=sparse(K(1:nnz_K,1),K(1:nnz_K,2),K(1:nnz_K,3),n_pd,n_pd);
    disp('M_sl is assembled');
    
% Build M_fo which is comprised by half discs and lateral chambers
M_fo = M_hd+M_lc;
disp('M_fo is assembled');

% Build M_ch which is the matrix to be used for integration
% over boundary domains occupied by channels.

%M_ch = M_fo + M_sl;
M_ch  = M_sl;
disp('M_ch is assembled');
    
%%%%%%%%%%%%%%%%%% Build Global Mass and Stiffness Matrices %%%%%%%%%%%%%%%
% Rewritten in terms of the mass and stiffness matrices this is
% cc_u*M_vol[u'] + 
% k_u*K_vol[u] +.5*beta_dark*nu_RE*epsilon_0*M_fo[u] =...
% - k*_{sigma,hyd}M_{sd}[E_st.*][u]+...
% .5*nu_RE*epsilon_0*M_fo[\alpha(v)]

% cc_v*M_vol[v'] + k_v*K_vol[v] =
% -1/(B_Ca*F)M_fo[J_ex(v)] + f_{Ca}/(2*B_Ca*F)M_fo[J_cG(u)]

% Together as one vector system, this reads
% [cc_u*M_vol zeros; zeros cc_v*M_vol][u';v'] +
% [k_u*K_vol + .5*beta_dark*nu_RE*epsilon_0*M_fo zeros;...
%  zeros k_v*K_vol][u;v] = 
% [-k*_{sigma,hyd}M_hd[E_st.*][u]+.5*nu_RE*epsilon_0*M_fo[\alpha(v)];...
%  -1/(B_Ca*F)M_fo[J_ex(v)] + f_{Ca}/(2*B_Ca*F)M_fo[J_cG(u)]]

% This conversion needed to pass volumic concentration
% of cGMP into the surface density. Rmbr E* concentration
% is only on hd.
vol2surfR= Volume_cone/Sigma_hd;

[I,J,V] = find(M_vol);
M_gl    = sparse([I;n_p3d+I], [J,n_p3d+J], [cc_u*V;cc_v*V],2*n_p3d,2*n_p3d);
disp('M_gl is assembled');

[I_Kvol,J_Kvol,V_Kvol] = find(K_vol);
[I_Mhd,J_Mhd,V_Mhd]    = find(M_hd);

%K_gl    = sparse([I_Kvol;I_Mfo;n_p3d+I_Kvol],[J_Kvol;J_Mfo;n_p3d+J_Mvol],...
%            [k_u*V_Kvol;.5*beta_dark*nu*epsilon_0*V_Mfo;k_v*V_Kvol],...
%            2*n_p3d,2*n_p3d);

K_gl    = sparse([I_Kvol;I_Mhd;n_p3d+I_Kvol],[J_Kvol;J_Mhd;n_p3d+J_Kvol],...
            [kk_u*V_Kvol;beta_dark*vol2surfR*V_Mhd;kk_v*V_Kvol],...
            2*n_p3d,2*n_p3d);
disp('K_gl is assembled');
    
save('Assembled_M_K','M_vol','K_vol','M_fo','M_hd','M_lc','M_pan','M_sl','M_ch',...
          'M_gl','K_gl','Sigma_cone','Sigma_sl','Sigma_hd','Volume_cone');
else
    load('Assembled_M_K');
    clear M_gl K_gl
    
    n_p3d = size(p_3d,2);
    vol2surfR= Volume_cone/Sigma_hd;

    [I,J,V] = find(M_vol);
    M_gl    = sparse([I;n_p3d+I], [J,n_p3d+J], [cc_u*V;cc_v*V],2*n_p3d,2*n_p3d);
    disp('M_gl is assembled');

    [I_Kvol,J_Kvol,V_Kvol] = find(K_vol);
    [I_Mhd,J_Mhd,V_Mhd]    = find(M_hd);
    
    K_gl    = sparse([I_Kvol;I_Mhd;n_p3d+I_Kvol],[J_Kvol;J_Mhd;n_p3d+J_Kvol],...
            [kk_u*V_Kvol;beta_dark*vol2surfR*V_Mhd;kk_v*V_Kvol],...
            2*n_p3d,2*n_p3d);
    disp('K_gl is assembled');
end

end

