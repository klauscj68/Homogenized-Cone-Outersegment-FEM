% costruisce la mesh
function [n_pd, n_int_pd, n_tri, n_sl, n_fo,  ...
         sl2pd, fo2pd, pd2int_pd, p_pd, t_pd, ...
         n_pd_sd, n_int_pd_sd, n_tri_sd, n_sl_sd, n_fo_sd, ...
         sl2pd_sd, fo2pd_sd, pd2int_pd_sd, p_pd_sd, t_pd_sd, ...
         Z_s, sd2sez, z_scaling]=...
     genemesh(R_b, R_t, H, theta_in, theta_fin, n_sez, flag_geom_sp, dz_0, ...
     n_sd, Z_sd, ...
     taglia, tol_R, tol_angle, n_ref_cyto, n_ref_id, ...
     plot_mesh, plot_num, inspect)
 
 

 
 
 
     
% output

% contatori
% n_pd è il numero di punti nel disco pivot
% n_int_pd è il numero di punti del disco pivot che non appartengono alla circonferenza esterna né ad incisure
% n_tri è il numero di triangoli in cui è partizionato il disco pivot
% n_os è il numero di punti sul perimetro del disco pivot
% n_p_inc è un vettore che contiene il numero di punti della traccia di ogni incisura sul disco pivot
% n_sez è il numero delle sezioni nel volume (assegnato in input)
% n_sd è il numero dei dischi speciali (assegnato in input)

% corrispondenze
% os2pd contiene gli indici dei nodi del bordo del disco pivot, in senso antiorario, nella numerazione del disco pivot
% inc2pd è un cell array: per ogni incisura contiene gli indici dei nodi, dal vertice al bordo, nella numerazione del disco pivot
% pd2int_pd associa ad ogni nodo sul disco pivot: 
%   zero, se il nodo è sulla buccia oppure su una incisura; 
%   il numero progressivo nella numerazione dei soli nodi del disco pivot non sulla buccia né su incisure, altrimenti
% sd2sez associa a ciascun disco speciale il numero di sezione cui esso è sovrapposto

% mesh
% p_pd è il vettore delle coordinate nodali del disco pivot
% t_pd è la matrice delle incidenze del disco pivot
% t_pd è la matrice delle incidenze del disco pivot inciso
% Z_s è il vettore delle quote delle sezioni rette

% a soli fini grafico-estetici in questa routine sono calcolati:
% p_vol, t_vol è la mesh di prismi retti a base triangolare nel volume
% p_sd, t_sd è la mesh di triangoli su uno qualsiasi dei dischi speciali
% p_os, t_os è la mesh di rettangoli sull'outer shell
% p_inc, t_inc sono strutture a cell array che contengono le mesh di rettangoli su ciascuna incisura
% nota: queste variabili non escono da genemesh, e non sono usate altrove

% nota bene:
% le variabili 
% n_pd, n_int_pd, n_tri, n_os, n_p_inc, ...
%          os2pd, inc2pd, pd2int_pd, p_pd, t_pd, ...
%          n_ipd, p_ipd, p_ipd_bea, t_ipd, ipd2pd
% escono anche nella forma 
%          n_pd_id, n_int_pd_id, n_tri_id, n_os_id, n_p_inc_id, ...
%          os2pd_id, inc2pd_id, pd2int_pd_id, p_pd_id, t_pd_id, ...
%          n_ipd_id, p_ipd_id, p_ipd_bea_id, t_ipd_id, ipd2pd_id
% che si riferisce alla mesh del disco con grado di raffittimento n_ref_id
% anziché n_ref_cyto
% le variabili senza suffisso _id servono per il problema della diffusione
% nel cytosol; le variabili col suffisso _id servono per il problema 
% della cascata sul disco


% numero di segmenti in cui è partizionato [0,H]

n_seg=n_sez-1;




% mesh del disco pivot, con e senza incisure
% nota bene: qui si passa n_ref_cyto, perché questo output serve per la mesh
% della diffusione hom
fprintf('\nGenerate the mesh for the special disc in the cytosol problem (mesh refinement %4i)\n',n_ref_cyto);
[n_pd, n_int_pd, n_tri, n_sl, n_fo, ...
        sl2pd, fo2pd, pd2int_pd, p_pd, t_pd]=...
    sezione_pivot(R_b, taglia, n_ref_cyto, tol_R, theta_in, theta_fin, tol_angle);



% informativa
fprintf('\nData relevant to the mesh for homogenized diffusion');
fprintf('\nmesh of the pivot disc %4i nodes; %4i triangles', n_pd, n_tri);
fprintf('\nMesh of the volume:                %4i nodes; %4i prismsi', n_pd*n_sez, n_tri*n_seg);
fprintf('\nMesh of the sliver:           %4i nodes; %4i rectangles', n_sl*n_sez, n_sl*n_seg);


% % stampe di controllo
% n_pd
% n_int_pd
% n_tri
% n_os
% os2pd
% pd2int_pd
% p_pd
% t_pd, 
% n_ipd
% p_ipd
% p_ipd_bea
% t_ipd
% ipd2pd
% for k=1:n_inc
%     'incisura', k
%     'n_p_inc(k)', n_p_inc(k)
%     'inc2pd{k}', inc2pd{k}
% end

% crea la mesh del disco inciso con grado di raffittimento n_ref_id
% infatti il numero di raffittimenti del problema della diffusione di E^* 
% sul disco inciso può essere maggiore di quello del problema omogeneizzato
fprintf('\nGenerate the mesh for the special disc in the disc problem (mesh refinement %4i)\n',n_ref_id);
[n_pd_sd, n_int_pd_sd, n_tri_sd, n_sl_sd, n_fo_sd,...
        sl2pd_sd, fo2pd_sd,  pd2int_pd_sd, p_pd_sd, t_pd_sd]=...
    sezione_pivot(R_b, taglia, n_ref_id, tol_R, theta_in, theta_fin, tol_angle);


% informativa
fprintf('\nMesh in the special disc: %i nodes; %i trianglesi\n', n_pd_sd, n_tri_sd);

% calcola le quote delle sezioni
fprintf('\nHeight of the cross sections\n');
[Z_s, sd2sez]=quote(H, n_sez, flag_geom_sp, dz_0, n_sd, Z_sd);



% Compute the scaling factor lambda to transform the cylinder in a cone

% tangent of the angle of the cone opening
tanalpha=(R_b-R_t)/H;

% Rescaling factor depending on the level z
z_scaling=1-Z_s*tanalpha/R_b;







if plot_mesh
	figure(1)
    hold on
% Mesh di triangoli nello special disk d
		l=patch('Vertices',p_pd','Faces',t_pd');    
        set(l,'facecolor','g')
        plot(p_pd(1,sl2pd),p_pd(2,sl2pd),'b','Linewidth',2.5)
    
    if plot_num
		for cont=1:n_pd
            aa=text(p_pd(1,cont),p_pd(2,cont),num2str(cont));
            set(aa,'color','b');
		end
		for cont=1:n_tri
            bb=text(mean(p_pd(1,t_pd(1:3,cont))),mean(p_pd(2,t_pd(1:3,cont))),num2str(cont));
            set(bb,'color','k');
		end
    end
	xlabel('x [\mu m]')
	ylabel('y [\mu m]')
    axis([-R_b*1.1 R_b*1.1 -R_b*1.1 R_b*1.1])
	view(2)
end





if plot_mesh
	figure(2)
    hold on
		l=patch('Vertices',p_pd_sd','Faces',t_pd_sd');    
        set(l,'facecolor','g')
    
        plot(p_pd_sd(1,sl2pd_sd),p_pd_sd(2,sl2pd_sd),'b','Linewidth',2.5)
        plot(p_pd_sd(1,fo2pd_sd),p_pd_sd(2,fo2pd_sd),'r','Linewidth',2.5)
    
    if plot_num
		for cont=1:n_pd_sd
            aa=text(p_pd_sd(1,cont),p_pd_sd(2,cont),num2str(cont));
            set(aa,'color','b');
		end
		for cont=1:n_tri_sd
            bb=text(mean(p_pd_sd(1,t_pd_sd(1:3,cont))),mean(p_pd_sd(2,t_pd_sd(1:3,cont))),num2str(cont));
            set(bb,'color','k');
		end
    end
	xlabel('x [\mu m]')
	ylabel('y [\mu m]')
    axis([-R_b*1.1 R_b*1.1 -R_b*1.1 R_b*1.1])
	view(2)
end







% mesh of prisms in the volume
% it is obtained by cartesian product of the pivot disc and the z quotes

% nodal coordinates: x e y in p_pd are replicated in row, one block for
% each z section
% z in Z_s is replicated in row element by element n_pd times 
p_vol=[ reshape(p_pd(1,:)'*z_scaling,1,n_pd*n_sez);...
        reshape(p_pd(2,:)'*z_scaling,1,n_pd*n_sez);...
        reshape(ones(n_pd,1)*Z_s,1,n_pd*n_sez)];
    
    
    

% adjacent matrix in the volume
t_vol=zeros(6,0);
for k=1:n_sez-1
    t_vol=[t_vol, [(k-1)*n_pd+t_pd; k*n_pd+t_pd]];
end

% plot mesh of prisms
if plot_mesh
	figure(3)
	% Mesh of prisms in the cone
    % bases
	l=patch('Vertices',p_vol','Faces',[t_vol([1 2 3],:)'; t_vol([4 5 6],:)']);
    % lateral surface
	l=patch('Vertices',p_vol','Faces',[t_vol([1 2 5 4],:)'; t_vol([2 3 6 5],:)'; t_vol([3 1 4 6],:)']);
	set(l,'facecolor','c')
	axis equal
	xlabel('x [\mu m]')
	ylabel('y [\mu m]')
	zlabel('z [\mu m]')
    axis([-R_b*1.1 R_b*1.1 -R_b*1.1 R_b*1.1 0 H])
	view(3)
end



% mesh of rectangles in the sliver
% they are obtained as cartesian product between the sliver trace in the
% pivot disc and the z quote vector

% nodal coordinates: x,y in p_sl is obtained from sl2pd, and replicated in
% row one block for each cross section
% z in Z_s is resplicated in row element by element n_sl times
% constructs the mesh
% coordinates x,y,z



p_sl=[ reshape(p_pd(1,sl2pd)'*z_scaling,1,n_sl*n_sez);...
        reshape(p_pd(2,sl2pd)'*z_scaling,1,n_sl*n_sez);...
        reshape(ones(n_sl,1)*Z_s,1,n_sl*n_sez)];
% incidenze
t_sl=zeros(4,0);
for k=1:n_sez-1
    t_sl=[t_sl, [(k-1)*n_sl+(1:n_sl); (k-1)*n_sl+[2:n_sl, 1]; k*n_sl+(1:n_sl); k*n_sl+[2:n_sl, 1]] ];
end

% disegna la mesh dei settori sull'os
if plot_mesh
	figure(4)
	% Mesh di rettangoli nell'outer shell
	l=patch('Vertices',p_sl','Faces',t_sl([1 2 4 3],:)');
	set(l,'facecolor','g')
	axis equal
	xlabel('x [\mu m]')
	ylabel('y [\mu m]')
	zlabel('z [\mu m]')
    axis([-R_b*1.1 R_b*1.1 -R_b*1.1 R_b*1.1 0 H])
	view(3)
end



presskey(inspect);



return
