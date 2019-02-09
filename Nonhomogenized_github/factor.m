% assembla le matrici globali e predispone per l'integrazione nel tempo
function [M_vol, M_fo,M_hd, M_lc, M_pan,M_sl, M_ch,...
    M_gl, K_gl,Sigma_cone,Sigma_sl,Sigma_hd,Volume_cone,...
    LL, UU, RR, p]=...
    factor(p_3d,f_3d,pr_3d,...
    n_gdl,...
    t_fin, n_step_t, metodo_cyto, theta, ...
    epsilon_0, nu, sigma, ...
    cc_u, kk_u, cc_v, kk_v)
%R_b, R_t, H, z_scaling, sd2sez, n_sez, n_sd, ...
%n_pd, n_tri, n_sl, sl2pd, n_fo, fo2pd, ...
%p_pd, t_pd, Z_s, ...
%gdl_vol, gdl_sd, gdl_sl,...

% M_gl, K_gl, LL, UU, RR, p are used for time integration
% M_vol, M_sd, M_sl are used for computing loads

% construction of mass and stiffness matrices in each region: vol, sd, sl
tic
%[M_vol, K_vol, M_sd, K_sd, M_sl, K_sl, M_fo]=...
%    assembla(n_pd, n_tri, n_sl, n_fo, n_sez, n_sd, ...                      % Adapt the p_pd scaling method
%    p_pd, t_pd, Z_s, sl2pd, fo2pd, sd2sez, z_scaling);
[M_vol, ~, M_fo,M_hd, M_lc, M_pan,M_sl, M_ch,...
          M_gl, K_gl, Sigma_cone, Sigma_sl, Sigma_hd,Volume_cone] =...
                      assembla( p_3d,f_3d,pr_3d );
toc

% % % % controlla l'assemblaggio
% % % check_MK(R, H, n_sez, n_sd, ...
% % %     n_pd, n_tri, n_os, n_inc, n_p_inc, ...
% % %     os2pd, inc2pd, ...
% % %     p_pd, t_pd, Z_s, Z_sd, inc, ...
% % %     M_vol, K_vol, M_sd, K_sd, M_os, K_os, M_inc, K_inc);

% assembla le matrici delle rigidezze e delle masse per il problema omogeneizzato
% mettendo insieme i pezzi di vol, sd, os, inc
%tic
%[K_gl, M_gl]=assembla_gl(R_b, R_t, H, n_pd, n_sl, ...                       % May need to change actually
%    n_sez, n_sd, ...
%    n_gdl, gdl_vol, gdl_sd, gdl_sl, ...
%    epsilon_0, nu, sigma, ...
%    cc_u, kk_u, cc_v, kk_v, ...
%    M_vol, K_vol, M_sd, K_sd, M_sl, K_sl);
%toc

% recupera un po' di memoria: le K non servono più, viceversa le M servono per i carichi
clear K_sd K_sl 

% passo di integrazione temporale
t_step=t_fin/n_step_t;

% fattorizzazione
if metodo_cyto(1)==1

    % permuta se richiesto
    tic
    if (metodo_cyto(3)==1) && (metodo_cyto(2)~=0)
        fprintf('\nSymmetric approximate minimum degree permutation\n');
        p = symamd(M_gl+(theta*t_step)*K_gl);
    elseif (metodo_cyto(3)==2) && (metodo_cyto(2)~=0)
        fprintf('\nSymmetric reverse Cuthill-McKee permutation\n');
        p = symrcm(M_gl+(theta*t_step)*K_gl);
    else
        fprintf('\nNo permutation\n');
        p=1:2*n_gdl;
    end
    toc

    % fattorizza se richiesto
    LL=[];
    UU=[];
    RR=[];
    tic
    switch metodo_cyto(2)
        case 1
            % fattorizzazione LU
            fprintf('\nFattorizzazione LU\n');
            [LL,UU]=lu(M_gl(p,p)+(theta*t_step)*K_gl(p,p));
        case 2
            % fattorizzazione di Cholesky
            fprintf('\nFattorizzazione di Cholesky\n');
            RR=chol(M_gl(p,p)+(theta*t_step)*K_gl(p,p));
        otherwise
            % non fattorizza
            fprintf('\nNo factorization\n');
    end
    toc

else
    % metodo ode
    % matrici vuote per evitare errore
    LL=[];
    UU=[];
    RR=[];
    p=[];
end

return
