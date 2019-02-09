% integrazione nel tempo
function [time, sol, curr_tot]=integrazione_tempo(p_3d,f_3d,...
    R_b, R_t, H, theta_in,theta_fin, n_sez, n_sd, ...
    n_gdl, ...
    t_fin, n_step_t, downsample, ...
    metodo_cyto, solver, theta, alpha, tol_fix, norma_inf, flag_Ca_clamp, ...
    M_vol, M_hd, M_sl, M_fo, M_ch, M_gl, K_gl, LL, UU, RR, p, ...
    n_sample, ...
    epsilon_0, nu, ...
    beta_dark, PDE_s, k_st, ...
    alpha_max, alpha_min, m_cyc, k_cyc, ... flag_ch,
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
    j_ex_sat, K_ex, ...
    u_ss, v_ss, E_st,...
    Sigma_cone, Sigma_sl, Sigma_hd,Volume_cone,...
    dof_vol,dof_sd, dof_sl, dof_fo, dof_ch,...
    dof_hd, dof_lc, dof_pan)



% passo di integrazione temporale
t_step=t_fin/n_step_t;
% istanti di calcolo
time=(0:n_step_t)'*t_step;

% storia della soluzione del problema per u, v relativa al primo campione
sol=zeros(2*n_gdl,n_step_t+1);
% storia della corrente totale relativa a tutti i campioni
curr_tot=zeros(n_sample,n_step_t/downsample+1);

% dato iniziale: prima u_ss, poi v_ss, ciascuno ripetuto n_gdl volte, in colonna
%sol0=[u_ss*ones(n_gdl,1); v_ss*ones(n_gdl,1)];
sol0=[u_ss; v_ss];
sol0=repmat(sol0,1,n_sample);

% Integrazione nel tempo
if metodo_cyto(1)==1
    % Integrazione con il metodo theta
    fprintf('\nIntegra nel tempo - Metodo Theta\n')

    % alloca il dato iniziale del primo campione
    sol(:,1)=sol0(:,1);
    % corrente dark
    [curr_tot(:,1)]=corrente( ...
    n_gdl, dof_ch, Sigma_sl,...
    sol0, ...
    R_b, R_t, H, theta_in,theta_fin, ... flag_ch,
    j_cg_max, m_cg, K_cg, ...
    j_ex_sat, K_ex, M_sl, M_hd, M_ch, nu, epsilon_0,...
    n_sample);


    % solin è inizializzata al dato iniziale
	solin=sol0;
	solfin=sol0;

    % alloca lo spazio per E_st_th, che conterrà il valore-theta della E*
    %E_st_th=sparse(n_gdl, n_sd);                                            % Dropped middle third n_sample index. Made sparse%
    
    tic
    
    % ciclo principale
    passo_downsample=1;
    rem_downsample=0;
    for passo=1:n_step_t
        save advance_cyto passo
        % istante finale dello step temporale
        t=time(passo+1);
        
		% mostra lo stato di avanzamento
		fprintf('Diffusione nel cytosol: istante %6.4f di %6.4f\n',t,t_fin);

        % valore theta di E_st
        %%%%%%%%%%%%% New
        Itot = [];%%%
        Jtot = [];%%%
        Vtot = [];%%%
        %%%%%%%%%%%%%
        for d=1:n_sd %Was E_st_th(:,d)= matrix expr in find
            f=(theta+rem_downsample)/downsample;
            %%%%%%%%%%%%%%%%%%%%% %New
            [I_ram,~,V_ram] = find(E_st{d}(:,passo_downsample)*(1-f)+E_st{d}(:,passo_downsample+1)*f);
            J_ram = zeros(size(I_ram,1),1)+d;%%%   %All this to be put in dth column
            Itot = [Itot;I_ram];%                  %Note that even tho we put diff sd's 
            Jtot = [Jtot;J_ram];%                  %in diff columns, the rows are indexed 
            Vtot = [Vtot;V_ram];%                  %by the global degrees of freedom. If findactch maps 
            %%%%%%%%%%%%%%%%%%%%%                  %to common disc should aggregate
% % %             E_st_th(:,1,d)=E_st{d}(:,1,passo)*(1-theta)+E_st{d}(:,1,passo+1)*theta;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %New
        E_st_th = sparse(Itot,Jtot,Vtot,n_gdl,n_sd);%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear I_ram J_ram V_ram Itot Jtot Vtot
		% calcolo di u e v alla fine dello step di integrazione, in applicazione del theta-metodo
		
		% soliter sarà l'incognita a fine step
		% come valore di tentativo, prendo l'incognita ad inizio step
		soliter=solin;
		
		% iterazione sul punto fisso (soliter)
		iter=0;
        again=true(1,n_sample);
		while any(again),
			
            % stima dell'incognita all'istante t_theta
            solth=solin(:,again)*(1-theta)+soliter(:,again)*theta;
            
            % separa il vettore delle u dal vettore delle v
            u_th=solth(1:n_gdl,:);
            v_th=solth(n_gdl+(1:n_gdl),:);
		
            % assembla i termini noti
            [f_th]=forcing_cyto(p_3d,f_3d,...
                            Sigma_cone, Sigma_sl, Sigma_hd,Volume_cone,...
                            R_b,R_t, H, theta_in, theta_fin, n_sez, n_sd, ...
                            epsilon_0, nu, ...
                            beta_dark, k_st, PDE_s,...
                            alpha_max, alpha_min, m_cyc, k_cyc, ...flag_ch,
                            B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
                            j_ex_sat, K_ex, ...
                            u_th, v_th, E_st_th, ...
                            M_sl, M_hd,...
                            dof_vol,dof_sd,dof_sl,dof_fo,dof_ch,...
                            dof_hd,dof_lc,dof_pan);

            % theta-metodo
            b=M_gl*solin(:,again) + t_step*(f_th - (1-theta)*(K_gl*solin(:,again)));
            switch metodo_cyto(2)
                case 0
                    solfin(:,again) = (M_gl+(theta*t_step)*K_gl) \ b;
                case 1
                    solfin(p,again) = UU\(LL \ b(p,:));
                case 2
                    solfin(p,again) = RR\(RR' \ b(p,:));
            end

            % rilassamento
            solfin(:,again)=(1-alpha)*soliter(:,again)+alpha*solfin(:,again);

            % clamp del calcio
            if flag_Ca_clamp
                % riporta la concentrazione del calcio in tutti i nodi
                % al valore dello steady-state
                solfin(n_gdl+1:end,again)=v_ss*ones(n_gdl,length(find(again)));
            end
            
            % fine se raggiunta la convergenza: norma L^inf
            if norma_inf
                % norma L^inf
%                 again= (max(abs(solfin-soliter),[],1)./max(abs(solfin),[],1))>tol_fix;
                again(again)= (max(abs(solfin(:,again)-soliter(:,again)),[],1)./max(abs(solfin(:,again)),[],1))>tol_fix;
            else
                % norma L^1
%                 again= (sum(abs(solfin-soliter),1)./sum(abs(solfin),1))>tol_fix;
                again(again)= (sum(abs(solfin(:,again)-soliter(:,again)),1)./sum(abs(solfin(:,again)),1))>tol_fix;
            end
            
            % per il passo successivo
            soliter=solfin;
            iter=iter+1;
    		fprintf('        iterazione %2i, campioni residui %4i\n',iter,length(find(again)));
            
		end % while punto fisso
        
        % per il passo temporale successivo
        solin=solfin;
    
        % salva i risultati sottocampionando
        rem_downsample=rem_downsample+1;
        if rem_downsample==downsample
            % aggiorna passo_downsample
            passo_downsample=passo_downsample+1;
            rem_downsample=0;

            % salva la storia della soluzione relativa al primo capione
            sol(:,passo_downsample)=solfin(:,1);
            
            % corrente totale
            [curr_tot(:,passo_downsample)]=corrente( ...
                            n_gdl, dof_ch, Sigma_sl,...
                            solfin, ...
                            R_b, R_t, H, theta_in,theta_fin, ...flag_ch,
                            j_cg_max, m_cg, K_cg, ...
                            j_ex_sat, K_ex, M_sl, M_hd, M_ch, nu, epsilon_0,...
                            n_sample);

        end

    end;
	toc

else
    error('ODE solver should not be method for time integration');
    % metodo ode
    %fprintf('\nIntegra nel tempo - Metodo %s\n',solver)
    %
    %tic
    %
    % matrice delle masse
    %options = odeset('Mass',M_gl,'InitialStep',t_step,'RelTol',tol_fix);
    %
    % ciclo su tutti i campioni
    %for samp=1:n_sample
    %    
    %    fprintf('Campione aleatorio %4i di %4i\n',samp,n_sample);
    %    
    %    % estrae E_st
    %    E_st_sample=cell(1,n_sd);
    %    for d=1:n_sd
    %        E_st_sample{d}=squeeze(E_st{d}(:,samp,:))';
    %    end
    %    
    %    % integra
    %    time_downsample=time(1:downsample:n_step_t+1);
    %    [~,sol_sample] = feval(solver, @odefun, time_downsample, sol0(:,samp), options,...
    %            R, H, n_sez, n_sd, ...
    %            n_pd, n_os, ...
    %            n_gdl, gdl_vol, gdl_sd, gdl_sl,...
    %            epsilon_0, nu, ...
    %            k_hyd, PDE_s, k_st, ...
    %            alpha_max, alpha_min, m_cyc, k_cyc, ...
    %            B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
    %            j_ex_sat, K_ex, ...
    %            K_gl, E_st_sample, t_step, n_step_t, t_fin, downsample, flag_Ca_clamp, ...
    %            M_vol, M_sd, M_os);
    %
    %    % salva la storia della soluzione relativa al primo capione
    %    if samp==1
    %        sol=sol_sample';
    %    end
    %
    %    % corrente totale
    %    [curr_tot(samp,:)]=corrente( ...
    %        n_p3d, dof_ch, Sigma_cone,...
    %        sol_sample, ...
    %        R_b, R_t, H, theta_in,theta_fin, ...
    %        flag_ch,j_cg_max, m_cg, K_cg, ...
    %        j_ex_sat, K_ex, M_sl, M_fo, M_ch, nu, epsilon_0);
    %    
    %end
end    
	toc

return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% usata da ode
% nell'integrare il sistema M_gl\dot{x}+K_gl{x}=f_th, ritorna f_th-K_gl{x}
%function y=odefun(t,x,...
%                R, H, n_sez, n_sd, ...
%                n_pd, n_os, ...
%                n_gdl, gdl_vol, gdl_sd, gdl_sl, ...
%                epsilon_0, nu, ...
%                k_hyd, PDE_s, k_st, ...
%                alpha_max, alpha_min, m_cyc, k_cyc, ...
%                B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
%                j_ex_sat, K_ex, ...
%                K_gl, E_st_sample, t_step, n_step_t, t_fin, downsample, flag_Ca_clamp, ...
%                M_vol, M_sd, M_os)
%
% mostra lo stato di avanzamento
%fprintf('Diffusione nel cytosol: istante %6.4f di %6.4f\n',t,t_fin);
%
% calcola E_st_sample all'istante t richiesto, con interpolazione lineare
%t_step=t_step*downsample;
%n_step_t=n_step_t/downsample;
%ind_t=fix(t/t_step)+1;
%ddt=rem(t,t_step);
%
% valore di E_st_sample all'istante t richiesto
%E_st_th=zeros(n_pd,n_sd);
%for d=1:n_sd,
%    if ind_t<n_step_t
%        E_st_th(:,d)=(E_st_sample{d}(ind_t,:)*(1-ddt/t_step)+E_st_sample{d}(ind_t+1,:)*ddt/t_step)';
%    else
%        E_st_th(:,d)=(E_st_sample{d}(n_step_t,:))';
%    end
%end
%
% estrae il vettore delle u dal vettore delle v
%u_th=x(1:n_gdl);
%v_th=x(n_gdl+(1:n_gdl));
%	
% assembla i termini noti
%[f_th]=carichi(R, H, n_sez, n_sd, ...
%    n_pd, n_os, ...
%    n_gdl, gdl_vol, gdl_sd, gdl_sl, ...
%    epsilon_0, nu, ...
%    k_hyd, PDE_s, k_st, ...
%    alpha_max, alpha_min, m_cyc, k_cyc, ...
%    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
%    j_ex_sat, K_ex, ...
%    u_th, v_th, E_st_th, ...
%    M_vol, M_sd, M_os);
%
% calcola il secondo membro del sistema
%y=-K_gl*x+f_th;
%
%if flag_Ca_clamp
    % clamp del calcio: annulla le derivate dei valori nodali di calcio
%    y(n_gdl+(1:n_gdl))=0;
%end
%return
