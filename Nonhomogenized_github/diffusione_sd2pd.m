% mappa la soluzione del problema per E_st_id nei dischi speciali incisi sui dischi non incisi
function [E_st_sample]=diffusione_sd2pd(n_pd, p_pd, p_pd_sd, t_pd_sd, ...
            n_sd, n_step_t, R_b, E_st_sd);
        
% I changed Giovanni's original script so that the E_st_sample now puts E*
% on dof_pvol. p_3d is new required input.
fprintf('\nProietta la diesterasi sulla mesh dei dischi speciali\n');



% proietta
    % la mesh utilizzata per la soluzione del problema di diffusione sul disco inciso (generata con n_ref_id raffittimenti)
    % è più fine di quella usata per il problema omogeneizzato (generata con n_ref_cyto raffittimenti), sicché
	% E_st_id{d}, soluzione del problema di diffusione sul disco inciso, 
	% deve essere proiettata sulla mesh del pivot disc della diffusione volumica, fornendo E_st_sample{d}
	
	% inizializza E_st_sample(d), con le stesse righe di E_st_sd{d}, ma con n_pd colonne
    E_st_sample=cell(1,n_sd);
    for d=1:n_sd
        E_st_sample(d)={zeros(size(E_st_sd{d},1),n_pd)};
    end
	
    % il fine è interpolare la soluzione sulla mesh rifinita per calcolarla
    % in ogni nodo p_pd della mesh di arrivo
    % per questo, trova il triangolo della mesh rifinita cui appartiene ogni nodo
    % della mesh di arrivo, quindi combina linearmente la soluzione della
    % mesh rifinita con le funzioni di forma campionate nel nodo della mesh
    % di arrivo
    % in pratica, avendo già cucito E_st_id, usa la mesh del pivot disk non
    % inciso, ma con raffittimento n_ref_id
    % in altri termini, usa p_pd_id, t_pd_id invece di p_ipd_id, t_ipd_id
    % tri e F hanno n_pd colonne
    [tri,F]=find_tri(p_pd_sd,t_pd_sd,p_pd,R_b);
    
    % combinazione lineare della soluzione sulla mesh raffittita con le funzioni F_elem
    % lista dei nodi dei triangoli della mesh raffittita che contengono i
    % punti p_pd, in colonna
    nodi=t_pd_sd(:,tri);
    nodi=nodi(1:3,:);                                                       %Colin added to cut out the 4th row which said all these triangles are in disc.
    nodi=nodi(:);
    % numero di campioni temporali
    % interpola linearmente i tre valori nodali di E_st_id{d}, per tutti i tempi
    for d=1:n_sd
        % soluzione calcolata nei nodi dei triangoli che contengono i punti p_pd, per tutti i tempi
        sol=E_st_sd{d}(:,nodi);                                             %Should be values on sd nodes
        % moltiplica la soluzione sol' per le repliche nei tempi del vettore
        % F, organizzato in colonna
        % questo fornisce una matrice che ha sulle righe i nodi a tre a
        % tre, sulle colonne i tempi
        % questa matrice è riordinata in modo da avere solo 3 righe
        % la somma sulle 3 righe esegue la combinazione lineare dei 3
        % valori nodali nei 3 verici del triangolo
        % il vettore riga risultante è poi riorganizzato, a gruppi di n_pd
        % elementi per ogni istante temporale
        E_st_sample{d}=reshape(sum(reshape(repmat(F(:),1,n_step_t+1).*sol',3,n_pd*(n_step_t+1)),1),n_pd,n_step_t+1)';
    end


return
