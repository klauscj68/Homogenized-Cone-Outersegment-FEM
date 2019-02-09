function [ dof_vol,dof_sd, dof_sl, dof_fo, dof_ch,...
           dof_hd, dof_lc, dof_pan] = dof(p_3d,f_3d)
% For all region types, output the degrees of freedom indexed as p_3d.
%   The regions are volume, special discs, sliver, folds,
%   horizontal discs, lateral chamber sides and panels. Outputs
%   are column vectors so easiest use with sparse
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
           j_ex_sat,K_ex,...
           flag_Ca_clamp] = data;

% dof_vol
n_p3d = size(p_3d,2);
dof_vol = cumsum(ones(n_p3d,1));

% dof_sd
tri = sd2sez(p_3d,f_3d);
ram = [];
for i=1:n_sd
    n_tri = size(tri{i},2);
    ram = [ram;reshape(tri{i},3*n_tri,1)];
end
dof_sd = unique(ram);
clear ram n_tri

% dof_lc
ram       = (f_3d(4,:)~=0).*(f_3d(5,:)==1); %chamber to cytosol
ram       = ram + (f_3d(5,:)==-2);          %sliver to cytosol
found     = find(ram)';
n_fd      = size(found,1);
ram = [];
for i=1:n_fd
    ram = [ram; reshape(f_3d(1:4,found(i)),4,1)];
end
dof_lc = unique(ram);
clear ram found n_fd

% dof_pan
ram       = (f_3d(5,:)==3)+(f_3d(5,:)==4);
found     = find(ram)';
n_fd      = size(found,1);
ram = [];
for i=1:n_fd
    ram = [ram; reshape(f_3d(1:4,found(i)),4,1)];
end
dof_pan = unique(ram);
clear ram found n_fd

% dof_sl
ram       = (f_3d(5,:)==5);
found     = find(ram)';
n_fd      = size(found,1);
ram = [];
for i=1:n_fd
    ram = [ram; reshape(f_3d(1:4,found(i)),4,1)];
end
dof_sl = unique(ram);
clear ram found n_fd

% dof_hd
found = find(f_3d(5,:)==-1);
n_fd  = size(found,2);
ram = [];
for i=1:n_fd
    ram = [ram; reshape(f_3d(1:3,found(i)),3,1)];
end
dof_hd = unique(ram);
clear ram found n_fd

% dof_fo
dof_fo = [dof_hd;dof_lc];
dof_fo = unique(dof_fo);

% dof_ch
%dof_ch  = [dof_fo;dof_sl];
%dof_ch  = unique(dof_ch);
dof_ch   = dof_sl;
end



function [tri] = sd2sez(p_3d,f_3d)
   %Output the triangular faces belonging to each special disc.  Note that
   %tri is a 1 by n_sd cell whose entries are 3 by n_tridsic matrices
   %listing the triangular elements of this special disc with respect to
   %p_3d.
[ ~,~,~,~,~,epsilon_0,nu,~,~,~,n_sd,~,~,~,n_sez,~,...
           ~] = data;
tri = cell(1,n_sd);
[Z_sd,~]=find_actch;
tol_z = nu*epsilon_0/(4*n_sez);%Height of smallest interdiscal chamber
                               %chamber is nu*eps0/2 which is broken up
                               %into n_sez.  If take half of the former
                               %divided by latter, no disc faces should
                               %be that close to the true faces. 
for i=1:n_sd
   % Find faces in a disc and which are tolerance
   % close to Z_sd(i)
   ram = (f_3d(4,:) == 0);        %Horizontal faces
   ram = ram.*(f_3d(5,:)==-1);   %Belonging to discs
   ram = ram.*...               %Inside z-tolerance
       (abs(p_3d(3,f_3d(1,:))...
        -Z_sd(i))<= tol_z);
   faces = find(ram);
   clear ram
   tri{i} = f_3d(1:3,faces);
   clear faces
end

end



