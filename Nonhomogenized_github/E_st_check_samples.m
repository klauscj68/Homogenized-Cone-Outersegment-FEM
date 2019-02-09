function [ sample ] = E_st_check_samples(p_3d,f_3d,...
                                         E_st,rad,angle)
%Evaluate E_st at given sample points. rad is rand(1,n) in [0,1] while
%angle is 2*pi*rand(1,n).


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

% mesh construction

        %Rebuild the cytosol disc cross section mesh
        [ p_pd_cyto,e_pd_cyto,t_pd_cyto,~,...
          n_pd_cyto,n_ed_cyto,n_tri_cyto,~,...
          tri_disc_cyto,tri_sl_cyto,...
          ~,~,~,~,~,...
           vert_disc_cyto,vert_sl_cyto,...
           ~,~,~,~,~]...
           = genmesh_cross(R_b,R_t,H,theta_in,theta_fin,...
                           epsilon_0,nu,sigma,...
                           taglia,tol_R,tol_angle,n_ref_cyto);

        % Trim the mesh down so it is now just the disc part
        % without the sliver.
        ram   = ismember(e_pd_cyto(1,:),vert_disc_cyto).*... %Find edges with both
        ismember(e_pd_cyto(2,:),vert_disc_cyto);     %vertices in the 
                                           %disc.
        ed_disc_cyto=find(ram);  

        p_pd_cyto = p_pd_cyto(:,vert_disc_cyto);    %Vertices, edges, and
        e_pd_cyto = e_pd_cyto(:,ed_disc_cyto);      %triangle lists trimmed
        t_pd_cyto = t_pd_cyto(:,tri_disc_cyto);     %to disc mesh.   

        n_pd_cyto = size(vert_disc_cyto,1);              %Update numbers of 
        n_tri_cyto = size(t_pd_cyto,2);                  %nodes, triangles,                 
        n_ed_cyto  = size(ed_disc_cyto,2);               %and edges.
                                     
        reind = @(i) find(vert_disc_cyto == i);          %Reindex the lists to
        e_pd_cyto(1:2,:) = arrayfun(reind,e_pd_cyto(1:2,:));  %reflect the new disc
        t_pd_cyto(1:3,:) = arrayfun(reind,t_pd_cyto(1:3,:));                   

                       
%%%%%%%%%%%%%%%%%%%%%%%%%

p_pd = p_pd_cyto(1:2,:);
t_pd = t_pd_cyto(1:3,:);

[~,dof_sd] = sd2sez(p_3d,f_3d,p_pd);

sample = cell(1,n_sd);
for kt=1:n_sd
% Is this z-scaling on p_pd
lambda = @(z) R_b/R_t - (R_b-R_t)/(H*R_t).*z;
z_scaling = lambda(Z_sd(kt));
p = z_scaling*p_pd;
R = z_scaling*R_t;
% Scale the p_find
p_find = ...
    R*[rad.*arrayfun(@(x) cos(x),angle);...
       rad.*arrayfun(@(x) sin(x),angle)];

[tri, F, ~, ~]=find_tri(p, t_pd, p_find, R);

%Grab this special disc
dof_sd_current = dof_sd{kt}; %ordered like p_pd
E_st_sd = E_st{kt}(dof_sd_current,:); %orders E_st values with p_pd index

n_pfind = size(p_find,2);
E_st_val = [];

for j=1:n_pfind
    %triangle pt belongs to
    triangle = tri(j);
    
    %triangle nodes
    nodes = t_pd(:,triangle);
    
    %E_st_sd all times
    E_st_loc = E_st_sd(nodes,:);
    %F is the shape functions eval'd at this point
    %E_st_loc are coefficients on shape functions
    %whereas F(:,i) are shape functions of that 
    %triangle evaluated at the point.
    ram    = (F(:,j)')*E_st_loc;
    E_st_val= [E_st_val;ram];
end
   sample{kt} = E_st_val;
end

end

function [tri,dof_sd] =...
          sd2sez(p_3d,f_3d,p_cyto) %p_pd_cyto(1:2,:) is what need
   %Output the triangular faces belonging to each special disc.  Note that
   %tri is a 1 by n_sd cell whose entries are 3 by n_tridsic matrices
   %listing the triangular elements of this special disc with respect to
   %p_3d. Also output the dof_sd node indices, w.r.t p_3d but ordered like
   %p_cyto as. p_cyto has the radius Rt while, p_3d has the cone geometry
   %radius. Alleged that dof_sd(i) = k means p_3d(1:2,k) is p_cyto(:,i) up
   %to homothety rescaling.
n_pcyto = size(p_cyto,2);     %number dof in sd
[ ~,~,~,~,~,epsilon_0,nu,~,~,~,n_sd,~,~,~,n_sez,~,...
           ~] = data;
tri  = cell(1,n_sd);          %Store triangles
dof_sd = cell(1,n_sd);        %Store where nodes belong

[Z_sd,~]=find_actch;
tol_z = nu*epsilon_0/(4*n_sez);%Height of smallest interdiscal chamber
                               %chamber is nu*eps0/2 which is broken up
                               %into n_sez.  If take half of the former
                               %divided by latter, no disc faces should
                               %be that close to the true faces. 
                               
% Reorganize p_cyto into a KD tree w.r.t x-y
% coordinates.  When match the p_3d nodes to
% p_cyto, we'll cycle over p_3d candidate nodes
% and compare to p_cyto by scanning over their
% KD tree. We must also rescale p_3d to the 
% common domain of p_cyto, so that we need
% only build 1 KD tree.
% Read geometric parameters
[ R_b,R_t,H,~,...
    ~] = data;
lambda = R_b/R_t - (R_b-R_t)/(H*R_t)*Z_sd;  %Scaling factor
lambdainv = 1./lambda;
                                            
KD_pcyto = KDTreeSearcher(p_cyto');         %p_cyto lives in disc R_t
disp('KD Tree Assembled')
for i=1:n_sd
   % Find faces in a disc and which are tolerance
   % close to Z_sd(i)
   ram = f_3d(4,:) == 0;         %Horizontal faces
   ram = ram.*(f_3d(5,:)==-1);   %Belonging to discs
   ram = ram.*...                %Inside z-tolerance
       (abs(p_3d(3,f_3d(1,:))...
        -Z_sd(i))<= tol_z);
   faces = find(ram);
   clear ram
   tri{i} = f_3d(1:3,faces);
   clear faces
   % Find where the nodes belong in relation to 
   % the tri{d} list. 
   ram = unique(tri{i}(:)); %All nodes in sd indexed 
                            %like dof_vol list.
   cand = p_3d(1:2,ram);    %Need to find where each column
                            %p_cyto is closest to which column
                            % cand. column of cand matches
                            % ram(column) index.
                            
   cand = lambdainv(i)*cand;%Rescale to domain of KD tree
   volissd = knnsearch(KD_pcyto,cand'); %Match vol ind to closest
                                        %sd node index, meaning 
                                        %vol ram(i) is sd volissd(i), ie
                                        %Ind -------> dof_sd
                                        %|   volissd
                                        %|ram
                                        %v
                                        %dof_vol
                                        %So map should be ram\circ
                                        %volissd^(-1)
   sdtovol = @(k) ram(find(volissd == k));
   dof_sd{i} = arrayfun(sdtovol,cumsum(ones(n_pcyto,1)));
end
end

