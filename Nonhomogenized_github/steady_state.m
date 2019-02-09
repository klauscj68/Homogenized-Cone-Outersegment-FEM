% soluzione numerica del sistema non lineare che descrive la soluzione stazionaria
function [u_ss,v_ss]=steady_state(Sigma_cone, Sigma_sl, Sigma_hd,Volume_cone,... new input Sigma_hd
    nu, epsilon_0, ...
    beta_dark, PDE_s, alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, j_ex_sat, K_ex, ...
    u_tent, v_tent, tol_stat)
% Hydrolysis term in assembla line 38 doesn't match 
% K_gl line 296. Also, the .5*nu*epsilon_0 term needs
% to be swapped out.   Input used to have k_hyd instead
% of beta_dark.

% This will eventually be used in place
% of the .5*nu*epsilon_0 terms which 
% adjust the surface density to the
% volume density.
%Sigma_ch = Sigma_cone - Sigma_sl;
Sigma_ch = Sigma_sl;
vol2surfR= Volume_cone/Sigma_hd;                                             %Careful since vol2surfR is conversion for E* machinery, actually want Volume_cone/Sigma_hd%

fprintf('\nSoluzione dello steady-state\n');

% dato iniziale
x0(1)=u_tent;
x0(2)=v_tent;

% risolve
x = fsolve(@fun, x0, optimset('TolFun', tol_stat), epsilon_0, nu, ...
    beta_dark, PDE_s, alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, j_ex_sat, K_ex,...
    Sigma_ch,vol2surfR); 

% mette la soluzione in u_ss e v_ss
u_ss=x(1);
v_ss=x(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fun(x, epsilon_0, nu, ...
    beta_dark, PDE_s, alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, j_ex_sat, K_ex,...
    Sigma_chv, vol2surf)
% estrae le incognite u,v da x
u=x(1);
v=x(2);

% output
y=zeros(1,2);

% prima equazione:
% bilancio dei flussi del cGMP sulle facce dei dischi
% The conversion ratio vol2surf passes volumic concentrations
% of cGMP to the surface density.
% This hydrolysis term doesn't match assembla line 296 in K_gl
%y(1)=k_hyd*PDE_s*u-(alpha_min+(alpha_max-alpha_min)*k_cyc^m_cyc/(k_cyc^m_cyc+v^m_cyc))*...
%    vol2surf;% Use to be (1/2*nu*epsilon_0);

y(1)=beta_dark*vol2surf*u-(alpha_min+(alpha_max-alpha_min)*k_cyc^m_cyc/(k_cyc^m_cyc+v^m_cyc))*...
    vol2surf;% Use to be (1/2*nu*epsilon_0);


% seconda equazione:
% bilancio dei flussi del Ca2+ sulla membrana plasmatica
% Appears to match terms in forcing_cyto
y(2)=j_ex_sat/(Sigma_chv*B_ca*F)*v/(K_ex+v)-j_cg_max*f_ca/(Sigma_chv*B_ca*F*2)*u^m_cg/(K_cg^m_cg+u^m_cg);
