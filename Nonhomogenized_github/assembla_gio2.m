% costruzione delle matrici di massa e rigidezza di ciascuna regione: vol, sd, os, inc
function [M_vol, K_vol, M_sd, K_sd, M_sl, K_sl, M_fo]=...
    assembla(n_pd, n_tri, n_sl, n_fo, n_sez, n_sd, ...
    p_pd, t_pd, Z_s, sl2pd, fo2pd, sd2sez, z_scaling)

% output

% M_vol, M_sd and M_os, are the mass matrices in the volume, sliver and special disc
% note: these matrices are not multiplied by any coefficient
% they are needed both to create the global matrices, 
% and the loading vector,
% the loading vector in each region is obtained as the product between the mass matrix relevant to the region and the nodal vector of the load 

fprintf('\nAssemble mass and stiffness matrices of vol, sd, sl\n')

% number of element layers
n_seg=n_sez-1;
n_arc_sl=n_sl-1;
n_arc_fo=n_fo-1;

% volume
% initialization
M_vol=zeros(36*n_tri*n_seg,3);                                               % Is this the right number?
K_vol=zeros(36*n_tri*n_seg,3);
% number of nonzero elements
nnz_K=0;
nnz_M=0;
% prisms in the volume
for s=1:n_seg
    for t=1:n_tri
        % indices of the triangular node in the pd numbering and of the
        % segment in the z-direction mesh
        ind_tri=t_pd(:,t);                                                   %Switch to pr_3d(:,t) and for loop there t%
        ind_sez=s:s+1;                                                       %Would automatically be ind_vol too.%
        % volume dof, in the progressive numbering of the volume             % local dof in the prism element?  
        ind_vol=reshape(repmat((ind_sez-1)*n_pd,3,1)+repmat(ind_tri,1,2),6,1); % Is it impt dof are sorted?
        % nodal coordinates in the cylindrical domain (will be rescaled to   % Point is to match the ref prism nodes 
        % get the coordinates in the cone domain)                             % to right index in physical space: ind_vol
        %X=p_pd(1,ind_tri); Not needed because now nodes rescaled            %Only question is is sort important 
        %Y=p_pd(2,ind_tri); By Matlab!
        %Z=Z_s(ind_sez);
        %lambda=z_scaling(ind_sez);
        % element matrices
        [K_elem]=prism_K(X,Y,Z,lambda);                                     % Replace with new routine
        [M_elem]=prism_M(X,Y,Z,lambda);
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
end;
% create M_vol and K_vol by adding contributes relevant to the same node    %They just mean you only have 3 lists until sparse makes matrix%
M_vol=sparse(M_vol(1:nnz_M,1),M_vol(1:nnz_M,2),M_vol(1:nnz_M,3),n_sez*n_pd,n_sez*n_pd);%And adds together entries belonging to same node pair%
K_vol=sparse(K_vol(1:nnz_K,1),K_vol(1:nnz_K,2),K_vol(1:nnz_K,3),n_sez*n_pd,n_sez*n_pd);


% special discs
% initialization of cell arrays
M_sd=cell(1:n_sd);
K_sd=cell(1:n_sd);

for d=1:n_sd
    M=zeros(9*n_tri,3);
    K=zeros(9*n_tri,3);
    % number of nonzero elements
    nnz_K=0;
    nnz_M=0;
    % triangles in the generic special disc
    for t=1:n_tri
        % nodal indices of the triangle (in the ordering of the pd mesh)
        ind_tri=t_pd(:,t);
        % actual nodal coordinate (take into account the scaling factor
        % relevant to the d-th special disc)
        X=z_scaling(sd2sez(d))*p_pd(1,ind_tri);                              %May have to redo.  I wrote something like this in main_id%
        Y=z_scaling(sd2sez(d))*p_pd(2,ind_tri);
        % element matrices
        [K_elem]=tri_K(X,Y);
        [M_elem]=tri_M(X,Y);
        % finds nonzero elements of K
        [I,J,V]=find(K_elem);
        nuovi_nnz=length(I);
        % stacks in K
        K(nnz_K+1:nnz_K+nuovi_nnz,:)=[ind_tri(I), ind_tri(J), V];
        nnz_K=nnz_K+nuovi_nnz;
        % finds nonzero elements of M
        [I,J,V]=find(M_elem);
        nuovi_nnz=length(I);
        % stacks in M
        M(nnz_M+1:nnz_M+nuovi_nnz,:)=[ind_tri(I), ind_tri(J), V];
        nnz_M=nnz_M+nuovi_nnz;
    end
    % somma i chiummi
    M_sd{d}=sparse(M(1:nnz_M,1),M(1:nnz_M,2),M(1:nnz_M,3),n_pd,n_pd);
    K_sd{d}=sparse(K(1:nnz_K,1),K(1:nnz_K,2),K(1:nnz_K,3),n_pd,n_pd);
end

% sliver                                                                      % Only mass matrix needed at bd
% initialization
M_sl=zeros(16*n_seg*n_sl,3);
K_sl=zeros(16*n_seg*n_sl,3);
% number of nonzero elements
nnz_K=0;
nnz_M=0;
% rectangles in the sliver
for s=1:n_seg
    for t=1:n_arc_sl
        % indices of segments nodes in the sliver and of of segments in the
        % mesh of the z direction
        ind_rett=t:t+1;
        ind_sez=s:s+1;
        % gdl of nodes in the sliver, in the progressive numbering of the
        % sliver
        % the nodes on the sliver are numbered in counterclock wise sense,
        % then from the bottom to the top
        ind_sl=reshape(repmat((ind_sez-1)*n_sl,2,1)+repmat(ind_rett',1,2),4,1);
        % coordinates of the base nodes of the restangule in the pivot disc
        % domain
        X=[p_pd(1,sl2pd(ind_rett(1))) p_pd(1,sl2pd(ind_rett(2)))];
        Y=[p_pd(2,sl2pd(ind_rett(1))) p_pd(2,sl2pd(ind_rett(2)))];
        Z=Z_s(ind_sez);
        lambda=z_scaling(ind_sez);
        % matrici di elemento
        [K_elem]=rectangule_K_gen(X,Y,Z,lambda);
        [M_elem]=rectangule_M_gen(X,Y,Z,lambda);
        % smembra K
        [I,J,V]=find(K_elem);
        nuovi_nnz=length(I);
        % inchiumma
        K_sl(nnz_K+1:nnz_K+nuovi_nnz,:)=[ind_sl(I), ind_sl(J), V];
        nnz_K=nnz_K+nuovi_nnz;
        % smembra M
        [I,J,V]=find(M_elem);
        nuovi_nnz=length(I);
        % inchiumma
        M_sl(nnz_M+1:nnz_M+nuovi_nnz,:)=[ind_sl(I), ind_sl(J), V];
        nnz_M=nnz_M+nuovi_nnz;
    end
end
% somma i chiummi
M_sl=sparse(M_sl(1:nnz_M,1),M_sl(1:nnz_M,2),M_sl(1:nnz_M,3),n_sl*n_sez,n_sl*n_sez);
K_sl=sparse(K_sl(1:nnz_K,1),K_sl(1:nnz_K,2),K_sl(1:nnz_K,3),n_sl*n_sez,n_sl*n_sez);




% siver
% initialization
M_fo=zeros(16*n_seg*n_fo,3);
% number of nonzero elements
nnz_M=0;
% rectangles in remaining part of the lateral surface (folds)
for s=1:n_seg
    for t=1:n_arc_fo
        % indices of segments nodes in the sliver and of of segments in the
        % mesh of the z direction
        ind_rett=t:t+1;
        ind_sez=s:s+1;
        % gdl of nodes in the folds, in the progressive numbering of the
        % folds
        % the nodes on the folds are numbered in counterclock wise sense,
        % then from the bottom to the top
        ind_fo=reshape(repmat((ind_sez-1)*n_fo,2,1)+repmat(ind_rett',1,2),4,1);
        % coordinates of the base nodes of the restangule in the pivot disc
        % domain
        X=[p_pd(1,fo2pd(ind_rett(1))) p_pd(1,fo2pd(ind_rett(2)))];
        Y=[p_pd(2,fo2pd(ind_rett(1))) p_pd(2,fo2pd(ind_rett(2)))];
        Z=Z_s(ind_sez);
        lambda=z_scaling(ind_sez);
        % matrici di elemento
        [M_elem]=rectangule_M_gen(X,Y,Z,lambda);
        % smembra M
        [I,J,V]=find(M_elem);
        nuovi_nnz=length(I);
        % inchiumma
        M_fo(nnz_M+1:nnz_M+nuovi_nnz,:)=[ind_fo(I), ind_fo(J), V];
        nnz_M=nnz_M+nuovi_nnz;
    end
end
% somma i chiummi
M_fo=sparse(M_fo(1:nnz_M,1),M_fo(1:nnz_M,2),M_fo(1:nnz_M,3),n_fo*n_sez,n_fo*n_sez);




return
