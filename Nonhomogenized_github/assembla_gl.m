% assembla le matrici delle rigidezze e delle masse per il problema omogeneizzato
% mettendo insieme i pezzi di vol, sd, os, inc
function [K_gl, M_gl]=assembla_gl(R_b, R_t, H, n_pd, n_sl, ...
    n_sez, n_sd, ...
    n_gdl, gdl_vol, gdl_sd, gdl_sl, ...
    epsilon_0, nu, sigma, ...
    cc_u, kk_u, cc_v, kk_v, ...
    M_vol, K_vol, M_sd, K_sd, M_sl, K_sl)



fprintf('\nAssemble global mass and stiffness matrices\n')

% initialize the matrices in which all non zero components will be stacked
K_gl=sparse(0,3);
M_gl=sparse(0,3);

% volume
% vectors of the volume dof's, ordered by columns
gdl_vol_col=reshape(gdl_vol,n_pd*n_sez,1);
% assembling
theta_0=1/(1+nu);
% finds non zeros elements of K_vol                                                 %We group these different dim sited M and K's%           
[I,J,V] = find(K_vol);                                                              %together bc they aggregated act on c' and c resp%
% stacks in K_gl                                                                    %Grouped on c' or c is what dictates build of M and K%
K_gl=[K_gl; [      gdl_vol_col(I),      gdl_vol_col(J),(1-theta_0)*kk_u*V]];        %Volume diffusion term for cGMP 
K_gl=[K_gl; [n_gdl+gdl_vol_col(I),n_gdl+gdl_vol_col(J),(1-theta_0)*kk_v*V]];        %and Ca respectively, notice index shift by #cGMP dof.
% find nonzero elements in M_vol
[I,J,V] = find(M_vol);
% stacks in M_gl
M_gl=[M_gl; [      gdl_vol_col(I),      gdl_vol_col(J),(1-theta_0)*cc_u*V]];       %Time mass matrix term in vol for cGMP
M_gl=[M_gl; [n_gdl+gdl_vol_col(I),n_gdl+gdl_vol_col(J),(1-theta_0)*cc_v*V]];       %and Ca respectively.


% special discs 
for d=1:n_sd
	% assembles
    gdl_sd_col=gdl_sd(:,d);
	% finds nonzero values in K_sd
	[I,J,V] = find(K_sd{d});
	% stacks in K_gl
	K_gl=[K_gl; [      gdl_sd_col(I),      gdl_sd_col(J),nu*epsilon_0*kk_u*V]];     %Special disc diffusion term for cGMP
	K_gl=[K_gl; [n_gdl+gdl_sd_col(I),n_gdl+gdl_sd_col(J),nu*epsilon_0*kk_v*V]];     %and Ca
	% find nonzero values in M_sd
	[I,J,V] = find(M_sd{d});
	% stacks in M_gl
	M_gl=[M_gl; [      gdl_sd_col(I),      gdl_sd_col(J),nu*epsilon_0*cc_u*V]];
	M_gl=[M_gl; [n_gdl+gdl_sd_col(I),n_gdl+gdl_sd_col(J),nu*epsilon_0*cc_v*V]];
end

% sliver


tanalpha=(R_b-R_t)/H;
cosalpha=1/sqrt(1+tanalpha^2);
% vectors of dofs in the sliver, ordered by column
gdl_sl_col=reshape(gdl_sl,n_sl*n_sez,1);
% assembles
% find nonzero elements K_sl
[I,J,V] = find(K_sl);
% stacks in K_gl
K_gl=[K_gl; [      gdl_sl_col(I),      gdl_sl_col(J),cosalpha*sigma*epsilon_0*kk_u*V]];
K_gl=[K_gl; [n_gdl+gdl_sl_col(I),n_gdl+gdl_sl_col(J),cosalpha*sigma*epsilon_0*kk_v*V]];
% find nonzero elements in M_sl
[I,J,V] = find(M_sl);
% stacks in M_gl
M_gl=[M_gl; [      gdl_sl_col(I),      gdl_sl_col(J),sigma*epsilon_0*cc_u*V]];
M_gl=[M_gl; [n_gdl+gdl_sl_col(I),n_gdl+gdl_sl_col(J),sigma*epsilon_0*cc_v*V]];



% somma finale dei chiummi
K_gl=sparse(K_gl(:,1),K_gl(:,2),K_gl(:,3),2*n_gdl,2*n_gdl);
M_gl=sparse(M_gl(:,1),M_gl(:,2),M_gl(:,3),2*n_gdl,2*n_gdl);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % APPROCCIO ALTERNATIVO, CHE DA' OUT OF MEMORY PER PROBLEMI GRANDI
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % inizializza
% K_gl=sparse(2*n_gdl,2*n_gdl);
% M_gl=sparse(2*n_gdl,2*n_gdl);
% 
% % volume
% % vettori dei gdl del volume, ordinati per colonna
% gdl_vol_col=reshape(gdl_vol,n_pd*n_sez,1);
% % assembla
% theta_0=1/(1+nu);
% K_gl(      gdl_vol_col,      gdl_vol_col) = K_gl(      gdl_vol_col,      gdl_vol_col)+(1-theta_0)*kk_u*K_vol;
% M_gl(      gdl_vol_col,      gdl_vol_col) = M_gl(      gdl_vol_col,      gdl_vol_col)+(1-theta_0)*cc_u*M_vol;
% K_gl(n_gdl+gdl_vol_col,n_gdl+gdl_vol_col) = K_gl(n_gdl+gdl_vol_col,n_gdl+gdl_vol_col)+(1-theta_0)*kk_v*K_vol;
% M_gl(n_gdl+gdl_vol_col,n_gdl+gdl_vol_col) = M_gl(n_gdl+gdl_vol_col,n_gdl+gdl_vol_col)+(1-theta_0)*cc_v*M_vol;
% 
% % dischi speciali (le matrici K_sd e M_sd sono le stesse per ogni disco speciale)
% for d=1:n_sd
% 	% assembla
%     gdl_sd_col=gdl_sd(:,d);
%     K_gl(      gdl_sd(:,d),      gdl_sd(:,d)) = K_gl(      gdl_sd(:,d),      gdl_sd(:,d))+nu*epsilon_0*kk_u*K_sd;
% 	M_gl(      gdl_sd(:,d),      gdl_sd(:,d)) = M_gl(      gdl_sd(:,d),      gdl_sd(:,d))+nu*epsilon_0*cc_u*M_sd;
% 	K_gl(n_gdl+gdl_sd(:,d),n_gdl+gdl_sd(:,d)) = K_gl(n_gdl+gdl_sd(:,d),n_gdl+gdl_sd(:,d))+nu*epsilon_0*kk_v*K_sd;
% 	M_gl(n_gdl+gdl_sd(:,d),n_gdl+gdl_sd(:,d)) = M_gl(n_gdl+gdl_sd(:,d),n_gdl+gdl_sd(:,d))+nu*epsilon_0*cc_v*M_sd;
% end
% 
% % outer shell
% % vettori dei gdl del'os, ordinati per colonna
% gdl_os_col=reshape(gdl_os,n_os*n_sez,1);
% % assembla
% K_gl(      gdl_os_col,      gdl_os_col) = K_gl(      gdl_os_col,      gdl_os_col)+sigma*epsilon_0*kk_u*K_os;
% M_gl(      gdl_os_col,      gdl_os_col) = M_gl(      gdl_os_col,      gdl_os_col)+sigma*epsilon_0*cc_u*M_os;
% K_gl(n_gdl+gdl_os_col,n_gdl+gdl_os_col) = K_gl(n_gdl+gdl_os_col,n_gdl+gdl_os_col)+sigma*epsilon_0*kk_v*K_os;
% M_gl(n_gdl+gdl_os_col,n_gdl+gdl_os_col) = M_gl(n_gdl+gdl_os_col,n_gdl+gdl_os_col)+sigma*epsilon_0*cc_v*M_os;
% 
% % incisure
% for m=1:n_inc
% 	% vettori dei gdl del'incisura m, ordinati per colonna
% 	gdl_inc_col=reshape(gdl_inc{m},n_p_inc(m)*n_sez,1);
% 	% assembla (inc(3,m) è l'angolo di beanza)
%     K_gl(      gdl_inc_col,      gdl_inc_col) = K_gl(      gdl_inc_col,      gdl_inc_col)+inc(3,m)*kk_u*K_inc{m};
% 	M_gl(      gdl_inc_col,      gdl_inc_col) = M_gl(      gdl_inc_col,      gdl_inc_col)+inc(3,m)*cc_u*M_inc{m};
% 	K_gl(n_gdl+gdl_inc_col,n_gdl+gdl_inc_col) = K_gl(n_gdl+gdl_inc_col,n_gdl+gdl_inc_col)+inc(3,m)*kk_v*K_inc{m};
% 	M_gl(n_gdl+gdl_inc_col,n_gdl+gdl_inc_col) = M_gl(n_gdl+gdl_inc_col,n_gdl+gdl_inc_col)+inc(3,m)*cc_v*M_inc{m};
% end

return
