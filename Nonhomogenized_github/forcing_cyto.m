function [ f_load ] = forcing_cyto(p_3d,f_3d,... New input Sigma_hd
    Sigma_cone, Sigma_sl, Sigma_hd,Volume_cone,...
    R_b,R_t, H, theta_in, theta_fin, n_sez, n_sd, ...
                            epsilon_0, nu, ...
                            Beta_dark, k_st, PDE_s,...
                            alpha_max, alpha_min, m_cyc, k_cyc, ...flag_ch,
                            B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
                            j_ex_sat, K_ex, ...
                            u_th, v_th, E_st_th, ...
                            M_sl, M_hd,...
                            dof_vol,dof_sd,dof_sl,dof_cyc,dof_ch,...
                            dof_hd,dof_lc,dof_pan)
% CHECK IF CONSTANTS IN FORCING TERM DATA ARE CORRECT
% IN CONE THE RIGHT CONVERSION BTWN
% SURFACE DENSITY AND VOLUME DENSITY IS NOT .5*EPS_0*NU 
% IN CONES.  WILL NEED TO USE Sigma_fo BELOW TO ADJUST.

Sigma_ch = Sigma_sl;
vol2surfR= Volume_cone/Sigma_hd;
dof_cyc = dof_hd;

%Build load vector for system cGMP and Ca diffusion
%   This script encodes the PDE contribution of cyclases synthesis of cGMP
%   as well as Ca2+ influx efflux by channels.  These are the nonlinear
%   terms of the forcing vector.  For moment hydrolysis and cyclase is put
%   at all the folds but not the panels nor sliver.  Channels are put
%   everywhere but the panels.  Note dof_cyc = dof_fo

% In particular, bc in main_id we map E_st onto that disc's nodes (vol 
% index) and keep it 0 everywhere else, can we not use M_hd rather than 
% M_sd?  Should be ok because those nodal functions are the only ones to be 
% nonzero on those special discs.

% Again Giovanni said actually you only want nu NOT nu_RE

% Together as one vector system, we had
% [cc_u*M_vol zeros; zeros cc_v*M_vol][u';v'] +
% [k_u*K_vol + .5*beta_dark*nu_RE*epsilon_0*M_fo zeros;...
%  zeros k_v*K_vol][u;v] = 
%  [-k*_{sigma,hyd}M_hd[E_st.*][u]+.5*nu_RE*epsilon_0*M_fo[\alpha(v)];...
%  -1/(B_Ca*F)M_fo[J_ex(v)] + f_{Ca}/(2*B_Ca*F)M_fo[J_cG(u)]]
% Actually put channels at folds & sliver! Cyclase still on folds.

% The forcing terms are on the right hand side of the equals.

% E_st_th is a 1 by n_sd cell.  Each entry is n_pd by n_sample by n_step_t
% matrix.  Because we're not doing random samples here n_sample = 1 and we
% should squeeze it.

%%%%%% First build alpha(v),J_ex(v), J_cG(v) as sparse column vec
n_p3d  = size(p_3d,2);
% Cyclase is only on horizontal discs
n_pcyc  = size(dof_cyc,1);
s_alpha = (alpha_min+(alpha_max-alpha_min)*...
          k_cyc^m_cyc./(k_cyc^m_cyc+v_th(dof_cyc).^m_cyc));
% Use ones because the alpha should be column vector
alpha   = sparse(dof_cyc,ones(n_pcyc,1),s_alpha,n_p3d,1);                         %In scripts I may have misused n_pd here

% Channels are on sliver 
n_pch   = size(dof_ch,1);
s_Jex   = j_ex_sat/(Sigma_ch)...*B_ca*F*(1+nu)*epsilon_0)...
          *v_th(dof_ch)./(K_ex+v_th(dof_ch));     %Sigma_cone part?
s_JcG   = j_cg_max/(Sigma_ch)*...f_ca/(2*B_ca*F*(1+nu)*epsilon_0))*...
          u_th(dof_ch).^m_cg./(K_cg^m_cg+u_th(dof_ch).^m_cg);
J_ex     = sparse(dof_ch,ones(n_pch,1),s_Jex,n_p3d,1);
J_cG     = sparse(dof_ch,ones(n_pch,1),s_JcG,n_p3d,1);

%%%%%%% Build the [E_st.*] sparse matrix 
% Actually don't need to build, already has been
[I_Est,~,V_Est] = find(E_st_th);                                             % In different function convert E_st into sparse column
% Diagonal matrix will just multiply each row by that entry                     % matrix over all dof_vol
M_Est = sparse(I_Est,I_Est,V_Est,n_p3d,n_p3d);                               % Note if findactch lumped different Z_sd's to the same
                                                                             % disc, it'll mean I_Est has repeats because they'll belong
                                                                             % to same global degree of freedom.  In turn sparse will lump
                                                                             % the contribution as should.

%%%%%%% Build the final forcing load column vector
%f_load = [-k_st*M_hd*M_Est*u_th+.5*nu*epsilon_0*M_fo*alpha;...
%          -1/(B_ca*F)*(M_fo+M_sl)*J_ex + f_ca/(2*B_ca*F)*(M_fo+M_sl)*J_cG];

% The vol2surfR converts the volumic concentration of
% cyclase into the surface density. This is appropriate
% because assuming the cyclase stays uniformly, steadily
% distributed even though Ca2+ drops.  That is the conversion
% acts on the alpha.  Don't believe need vol2surfR acting on
% the k_st term because that conversion has only been used to 
% pass [PDE]_vol to [PDE]_surf in the dark activated case (of
% course also for cyclase) while M_Est is already the 
% diffusion on the surface disc.

%f_load = [-k_st*M_hd*M_Est*u_th+vol2surfR*M_fo*alpha;...
%          -1/(B_ca*F)*(M_fo+M_sl)*J_ex + f_ca/(2*B_ca*F)*(M_fo+M_sl)*J_cG];

% Should be ok to use M_fo instead of M_hd because just
% like with M_hd, point is that M_Est is only supp'd
% at the boundary disc faces. Only need M_sl in Ca part for channels
f_load = [-k_st*M_hd*M_Est*u_th+vol2surfR*M_hd*alpha;...
          -1/(B_ca*F)*(M_sl)*J_ex + f_ca/(2*B_ca*F)*(M_sl)*J_cG];
      
end



