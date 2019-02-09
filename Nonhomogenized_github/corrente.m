% corrente totale, per tutti i campioni
function [curr_tot]=corrente( ...
    n_p3d, dof_ch, Sigma_sl,...
    solfin, ...
    R_b, R_t, H, theta_in,theta_fin, ...flag_ch,
    j_cg_max, m_cg, K_cg, ...
    j_ex_sat, K_ex, M_sl, M_hd, M_ch, nu, epsilon_0,...
    n_sample)

% Need n_p3d and Sigma_cone.

% solfin contains solition at a fixed time
% each column refers to a different stocastich sample 
curr_tot=zeros(1,n_sample);
% i carichi sono interpolati in accordo alla rappresentazione tramite le funzioni di forma
%Sigma=zeros(1,3);
% Lateral surface of the sliver                                                %Need to change how areas handled%
%tanalpha=(R_b-R_t)/H;
%eta=tanalpha/R_b;
%Sigma(1)=(theta_fin-theta_in)*R_b*sqrt(1+tanalpha^2)*(H-1/2*eta*H^2);
%
% Total disc membrane surface
%Sigma(2) = faces_area(R_t,R_b,H,nu,epsilon_0);






% Total esternal surface of the folds
%Sigma(3)=(2*pi-theta_fin+theta_in)*R_b*sqrt(1+tanalpha^2)*(H-1/2*eta*H^2);


% Compute surface where channels are located (uniformly distributed)
%Sigma_cone=sum(flag_ch.*Sigma);



% In nonhomogenized case, our channels won't have been
% homogenized into the interior.
%if flag_ch(2)==true
%    
    % solution in the volume
    % gdl_vol ordered by column
%    gdl_vol_col=reshape(gdl_vol,n_pd*n_sez,1);
    % estrae la soluzione sui nodi di os
%    u_vol=solfin(      gdl_vol_col,:);
%    v_vol=solfin(n_gdl+gdl_vol_col,:);
%
    % storia dei flussi di corrente nei nodi
%    dens_curr_cGMP=2*j_cg_max/(Sigma_cone*(1+nu)*epsilon_0)*(u_vol.^m_cg)./(K_cg^m_cg+u_vol.^m_cg); % Do we need (1+nu)*epsilon_0?
%    dens_curr_Ca  =2*j_ex_sat/(Sigma_cone*(1+nu)*epsilon_0)*v_vol./(K_ex+v_vol);
%
    % integrates on the volume
%    curr_cGMP=(ones(1,n_pd*n_sez)*M_vol*dens_curr_cGMP);                    %Need to change this to where channels are, ie M_fo+M_sl-->M_vol%
%    curr_Ca  =(ones(1,n_pd*n_sez)*M_vol*dens_curr_Ca);
%
    % corrente totale, per tutti i campioni
%    curr_tot=curr_tot+curr_cGMP+curr_Ca;
%end



%if flag_ch(1)==true                                                         % Do this case but M_fo+M_sl%
    % storia della soluzione sull'outer shell
    % gdl_os ordinato in colonna
    % gdl_sl_col=reshape(gdl_sl,n_sl*n_sez,1);                                % If needed do dof_sl%
    % estrae la soluzione sui nodi di os
                                     %The 2nd index used to be n_sample
                                     %which is my code is 1.
    %u_ch=solfin(dof_ch);            %u_sl=solfin(gdl_sl_col,:);      
    %v_ch=solfin(n_p3d+dof_ch);      %Double check the indexing is right
    u_ch = sparse(dof_ch,ones(size(dof_ch,1),1),...
           solfin(dof_ch),n_p3d,1);
    v_ch = sparse(dof_ch,ones(size(dof_ch,1),1),...
           solfin(n_p3d+dof_ch),n_p3d,1);   
    
    
    % storia dei flussi di corrente nei nodi
    dens_curr_cGMP=j_cg_max/Sigma_sl*(u_ch.^m_cg)./(K_cg^m_cg+u_ch.^m_cg); % Rewrite as sparse matrix n_p3d size & nonzero at dof_sl 
    dens_curr_Ca  =j_ex_sat/Sigma_sl*v_ch./(K_ex+v_ch);                    % "

    % integra sulla superficie laterale
    curr_cGMP=(ones(1,n_p3d)*(M_ch*dens_curr_cGMP));                      % ones(1,n_p3d)*M_sl*dens_curr_cGMP%
    curr_Ca  =(ones(1,n_p3d)*(M_ch*dens_curr_Ca));                        % "%

    % corrente totale, per tutti i campioni
    curr_tot=curr_tot+curr_cGMP+curr_Ca;


%if flag_ch(3)==true
    % storia della soluzione sull'outer shell
    % gdl_os ordinato in colonna
%    gdl_fo_col=reshape(gdl_fo,n_fo*n_sez,1);
    % estrae la soluzione sui nodi di os
%    u_fo=solfin(      gdl_fo_col,:);
%    v_fo=solfin(n_gdl+gdl_fo_col,:);

    % storia dei flussi di corrente nei nodi
%    dens_curr_cGMP=j_cg_max/Sigma_cone*(u_fo.^m_cg)./(K_cg^m_cg+u_fo.^m_cg);
%    dens_curr_Ca  =j_ex_sat/Sigma_cone*v_fo./(K_ex+v_fo);

    % integra sulla superficie laterale
%    curr_cGMP=(ones(1,n_fo*n_sez)*M_fo*dens_curr_cGMP);
%    curr_Ca  =(ones(1,n_fo*n_sez)*M_fo*dens_curr_Ca);

    % corrente totale, per tutti i campioni
%    curr_tot=curr_tot+curr_cGMP+curr_Ca;
%end

return

end



