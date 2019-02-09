function [ p_3d,f_3d,pr_3d ] = genmesh
%Master script to generate the mesh both on discs and in incised volume.
%   This routine is called by main to generate the meshes both on the
%   activated discs as well as in the interior volume. It has an option
%   value "plotmesh' to indicate whether it's desired that the meshes be
%   plotted or not.  Note this routine builds most of the mesh on a 
%   cylinder and then at the last scales it to a cone. Also, when Giovanni 
%   Bisegna wrote their code they sorted the nodes, using 'sort', in a
%   counter-clockwise manner and increasing radius.

%   p_3d is a (5xn_p3d) list.  The first 3 rows are the x,y,z coordinates
%   of the nodes. The fourth row says whether it belongs to the disc or
%   sliver, and the 5th row indicates which cross section arc it may
%   belong. 

%   f_3d is a (5xn_f3d) list.  Its first four rows encode nodal vertices
%   comprising the face listed in a continuous tour (either clockwise or
%   counter-clockwise along the face).  If it is triangular rather than
%   rectangular its fourth vertex is shown as 0.  Its fifth row, when
%   pertaining to a triangular element, gives 1 or 2 to say whether in disc
%   or sliver.  If pertaining to a rectangular element, it gives a number 0
%   to 5 to say whether it belongs to an exterior lateral surface matched
%   to the arc edge 1 to 5.

% pr_3d is a (7xn_pr3d) list.  Its first 6 rows are the node indices
% comprising the prism.  Its last row encodes whether it belongs to the
% disc or the sliver by encoding a 1 or 2, like the triangles list.


%Make a call to geometric parameters defining the problem
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

%Create the 2d cross section mesh with genmesh_cross
[ p_pd,e_pd,t_pd,e_arcpd,n_pd,n_ed,n_tri,n_arced,...
           tri_disc,tri_sl,...
           edgarc1,edgarc2,edgarc3,edgarc4,edgarc5,...
           vert_disc,vert_sl,...
           vertarc1,vertarc2,vertarc3,vertarc4,vertarc5]...
           = genmesh_cross(R_b,R_t,H,theta_in,theta_fin,...
                           epsilon_0,nu,sigma,...
                           taglia,tol_R,tol_angle,n_ref_cyto);  

%%%%%%%%%%%%%%%%  Build 3D Mesh in Cytosol B4 Cone Scaling %%%%%%%%%%%%%%%%
%This part more or less accomplishes what quote did in Giovanni's code.
%The code is in blocks.  First block builds a single 3d interdiscal
%chamber. The 2nd block then stacks n_chambers many of these, and the third
%block connects the sliver region between chambers.

%%%%% Build a model interdiscal chamber of unit height by
%%%%% stacking n_sez many cross sections and then building
%%%%% prisms by connecting up and down. Observe that this 
%%%%% chamber will join the sliver connects at the 
%%%%% 'vert_sl' nodes and the '(n_sez-1)*n_pd + vert_sl'
%%%%%  nodes. The joined edge indices are 'e_sl' and 
%%%%% '(n_sez-1)*n_ed + e_sl'. Similarly, the 
%%%%% triangles are 'tri_sl' and
%%%%% '(n_sez-1)*n_tri + tri_sl'.
%%%%% Finally, in top and bottom faces of chamber which
%%%%% are also in the disc, modify their 1's to be 4 and
%%%%% 3 respectively

zaxis = 0:1/(n_sez-1):1;            %Chamber has n_sez many joined cross
                                    %sections. zaxis gives the heights at
                                    %which cross sections are to be placed
                                    
p_ch = [zeros(1,n_pd);p_pd];              %Prepare to duplicate the
                                          %vertices listed as rows 
                                          %lists. This means
                                          %encoding a z-coordinate which we
                                          %do in the TOP row because of
                                          %membership information given in
                                          %genmesh_cross below.
                                                  
e_ch = e_pd;                              %Prepare to duplicate the
                                          %edges indexed by
                                          %vertices and by edge
                                          %number lists.

t_ch = t_pd;                              %Prepare to duplicate the
                                          %triangles listed by
                                          %vertex number lists. But also
                                          %keep subdomain membership bc
                                          %will use to say whether prism is
                                          %in domain or sliver or which
                                          %cross section triangles admit
                                          %the signaling biochemistry.

%This duplicates all the lists for each
%new cross section without creating 
%anything new.
for i=1:n_sez-1        
    p_ch = [p_ch [zaxis(1,i+1)*...                       %Concatenate from 
                   [ones(1,n_pd);zeros(4,n_pd)]+...      %left to right the
                   [zeros(1,n_pd);p_pd]]];               %nodes at new z
                                                         %heights.                                                     
    e_ch = [e_ch ...
                     [i*n_pd+e_pd(1:2,:);...          %Concatenate from
                      e_pd(3,:)]...                   %left to right the 
            ];                                        %vertex indices their
                                                      %way and the edge
                                                      %indices theirs. The
                                                      %vertex indices match
                                                      %p_pd_chamber because
                                                      %that was tiled left
                                                      %to right in n_pd
                                                      %blocks.
                                                      
    t_ch = [t_ch ...
            [i*n_pd+t_pd(1:3,:);t_pd(4,:)]...
            ];                                        %Shift the vertex 
                                                      %indices which define
                                                      %the triangles while
                                                      %concatenating left
                                                      %to right.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This was added to attempt to identify where boundary faces are
% Modify the top and bottom tri's which are 
% in the disc to encode that they belong 
% there by writing -1 instead of 1.

bottom = find(t_ch(4,1:n_tri) == 1);
top    = (n_sez-1)*n_tri + find(t_ch(4,(n_sez-1)*n_tri+1:n_sez*n_tri)...
              == 1);
t_ch(4,[bottom top]) = -1;
clear bottom top
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Connect cross sections along z.  This involves 
% augmenting the current edges lists and building 
% a new faces and prisms list.

% Build a faces list by connecting edges that are
% one up and one down. The first n_sez*n_tri 
% columns are triangles on cross sections the 
% last n_sez*n_ed are the vertical faces.
f_ch =         [t_ch(1:3,:);zeros(1,n_sez*n_tri);...
                t_ch(4,:)];
                                                %The triangles are already 
                                                %faces, but they only have
                                                %three edges whereas the
                                                %vertical parts of prisms
                                                %will be defined by four.
                                                                           
                                                %When add in edge built
                                                %prisms, the 5th row
                                                %instead of interior/sliver
                                                %will encode boundary arc
                                                %relationship.  To know if
                                                %5th row index is interior/
                                                %sliver or a geometric arc,
                                                %check if 4th row is 0.
for i=1:n_sez-1
    Permutation = [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];
    Permutation = Permutation*...                  %Permute the nodes so 
                  [...                             %the face has an 
                   e_ch(1:2,(i-1)*n_ed+1:i*n_ed);  %oriented tour of nodes.
                   e_ch(1:2,i*n_ed+1:(i+1)*n_ed)...%This connects the ith 
                  ];                               %face to the (i+1)th 
                                                   %face.
    f_ch = [f_ch ...                                    %Augment the tri's
                    [Permutation;...                    %with the vertices
                     e_pd(3,:)]...                      %that belong to the
           ];                                           %top and bottom 
                                                        %edges. Also, order
                                                        %nodes on face so
                                                        %if map std square
                                                        %to face in
                                                        %clockwise or
                                                        %counter-clockwise
                                                        %manner there is no
                                                        %twist, ie so the
                                                        %nodes on face are
                                                        %(counter)clockwise
                                                        %The zeros record
                                                        %unknown volume
                                                        %region membership
                                                        %while last 5 rows
                                                        %of e_ch record arc
                                                        %membership. Each
                                                        %block underneath
                                                        %Permutation should
                                                        %just be the
                                                        %repeating edge arc
                                                        %membership within
                                                        %each cross
                                                        %section. Hence
                                                        %n_ed sizes.                                                                                                       
end
%Build prism list by connecting triangles that are
%one up and one down.
pr_ch = [];                                 %Each column will be a prism 
                                            %that is itself defined by six
                                            %vertices.
for i=1:n_sez-1
    pr_ch = [pr_ch ...                                           %Connect 
                      [t_ch(1:3,(i-1)*n_tri+1:i*n_tri);...       %triangles
                       t_ch(:,i*n_tri+1:(i+1)*n_tri)]...         %that are 
            ];                                                   %directly
                                                      %below and above
                                                      %which lie in a
                                                      %common volume
                                                      %region, hence only
                                                      %need the second call
                                                      %of t_ch's membership
                                                      %rows. Vertices are
                                                      %ordered
                                                      %counter-clockwise
                                                      %bottom then
                                                      %directly above's, eg
                                                      %a,b,c,a',b',c'.
                                                               
end
                                                


% Augment edges list by connecting vertices that
% are difference height 1 in up and down faces.
for i=1:n_sez-1
    e_ch = [e_ch ...
                     [(i-1)*n_ed+1:i*n_ed;...       %lower face nodes
                      i*n_ed+1:(i+1)*n_ed;...       %upper face nodes
                      e_pd(3,:)]...                 %store the arcs they
                     ];                             %may belong to.
                                                                       
end
n_pch  = size(p_ch,2);                  %Update number of nodes, edges
n_ech  = size(e_ch,2);                  %faces and prisms.
n_fch  = size(f_ch,2);
n_prch = size(pr_ch,2); 

%%%%% Build a model unit sliver chamber
p_sl = p_pd(:,vert_sl);                 %Record vertices in sliver
n_psl= size(vert_sl,1);

ram = ismember(e_pd(1:2,:),vert_sl);    %Find edges with both vertices
ram = ram(1,:).*ram(2,:);               %belonging to the sliver.
ram = find(ram);                        %Grab columns with such edges

e_sl = e_pd(:,ram);                     %Sliver edges given by old vertex
n_esl= size(e_sl,2);                    %indexing system of p_pd.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This was added to indentify the exterior sliver chamber's patch boundary
% when coming to assemble mass and stiffness matrices. This gets it copied
% into the f_slch when time comes or should.
ram = find(e_sl(3,:)==2);
e_sl(3,ram) = -2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear ram
t_sl = t_pd(:,tri_sl);                  %Grab triangles in sliver still 
                                        %indexed by old vertex system.
ram_e = e_sl;
ram_t = t_sl;

for i=1:n_psl
    ev = find([e_sl(1:2,:);...
           zeros(1,size(e_sl,2))] == vert_sl(i,1));
    tv = find([t_sl(1:3,:);...
           zeros(1,size(t_sl,2))] == vert_sl(i,1));
    ram_e(ev) = i;
    ram_t(tv) = i;
end

e_sl = ram_e;                            %e_sl and t_sl are now reindexed
t_sl = ram_t;                            %to match p_sl
clear ram_e ram_t
                     
n_tsl=size(t_sl,2);

                                    
p_slch = [zeros(1,n_psl);p_sl];           %Prepare to duplicate the
                                          %vertices listed as rows 
                                          %lists. This means
                                          %encoding a z-coordinate which we
                                          %do in the TOP row because of
                                          %membership information given in
                                          %genmesh_cross below.
                                                  
e_slch = e_sl;                            %Prepare to duplicate the
                                          %edges indexed by
                                          %vertices and by edge
                                          %number lists.

t_slch = t_sl;                            %Prepare to duplicate the
                                          %triangles listed by
                                          %vertex number lists. But also
                                          %keep subdomain membership bc
                                          %will use to say whether prism is
                                          %in domain or sliver or which
                                          %cross section triangles admit
                                          %the signaling biochemistry.

%This duplicates all the lists for each
%new cross section without creating 
%anything new.
for i=1:n_sez-1        
    p_slch = [p_slch [zaxis(1,i+1)*...                   %Concatenate from 
                   [ones(1,n_psl);zeros(4,n_psl)]+...    %left to right the
                   [zeros(1,n_psl);p_sl]]];             %nodes at new z
                                                         %heights.                                                     
    e_slch = [e_slch ...
                     [i*n_psl+e_sl(1:2,:);...         %Concatenate from
                      e_sl(3,:)]...                   %left to right the 
            ];                                        %vertex indices their
                                                      %way and the edge
                                                      %indices theirs. The
                                                      %vertex indices match
                                                      %p_pd_chamber because
                                                      %that was tiled left
                                                      %to right in n_pd
                                                      %blocks.
                                                      
    t_slch = [t_slch ...
            [i*n_psl+t_sl(1:3,:);t_sl(4,:)]...
            ];                                        %Shift the vertex 
                                                      %indices which define
                                                      %the triangles while
                                                      %concatenating left
                                                      %to right.
end

% Connect cross sections along z.  This involves 
% augmenting the current edges lists and building 
% a new faces and prisms list.

% Build a faces list by connecting edges that are
% one up and one down. The first n_sez*n_tri 
% columns are triangles on cross sections the 
% last n_sez*n_ed are the vertical faces.
f_slch =         [t_slch(1:3,:);zeros(1,n_sez*n_tsl);...
                  t_slch(4,:)];
                                                %The triangles are already 
                                                %faces, but they only have
                                                %three edges whereas the
                                                %vertical parts of prisms
                                                %will be defined by four.
                                                %Last row is to store
                                                %membership of the top or
                                                %bottom face of the
                                                %chamber. 1 is bottom, 2 is
                                                %top.
                                                %When add in edge built
                                                %prisms, the 5th row
                                                %instead of interior/sliver
                                                %will encode boundary arc
                                                %relationship.  To know if
                                                %5th row index is interior/
                                                %sliver or a geometric arc,
                                                %check if 4th row is 0.
for i=1:n_sez-1
    Permutation = [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];
    Permutation = Permutation*...                  %Permute the nodes so 
                  [...                             %the face has an 
                   e_slch(1:2,(i-1)*n_esl+1:i*n_esl);  %oriented tour of nodes.
                   e_slch(1:2,i*n_esl+1:(i+1)*n_esl)...%This connects the ith 
                  ];                               %face to the (i+1)th 
                                                   %face.
    f_slch = [f_slch ...                                %Augment the tri's
                    [Permutation;...                    %with the vertices
                     e_sl(3,:)]...                      %that belong to the
           ];                                           %top and bottom 
                                                        %edges. Also, order
                                                        %nodes on face so
                                                        %if map std square
                                                        %to face in
                                                        %clockwise or
                                                        %counter-clockwise
                                                        %manner there is no
                                                        %twist, ie so the
                                                        %nodes on face are
                                                        %(counter)clockwise
                                                        %The zeros record
                                                        %unknown volume
                                                        %region membership
                                                        %while last 5 rows
                                                        %of e_ch record arc
                                                        %membership. Each
                                                        %block underneath
                                                        %Permutation should
                                                        %just be the
                                                        %repeating edge arc
                                                        %membership within
                                                        %each cross
                                                        %section. Hence
                                                        %n_ed sizes.                                                                                                       
end
%Build prism list by connecting triangles that are
%one up and one down.
pr_slch = [];                               %Each column will be a prism 
                                            %that is itself defined by six
                                            %vertices.
for i=1:n_sez-1
    pr_slch = [pr_slch ...                                   %Connect 
                      [t_slch(1:3,(i-1)*n_tsl+1:i*n_tsl);... %triangles
                       t_slch(:,i*n_tsl+1:(i+1)*n_tsl)]...   %that are 
            ];                                               %directly
                                                      %below and above
                                                      %which lie in a
                                                      %common volume
                                                      %region, hence only
                                                      %need the second call
                                                      %of t_ch's membership
                                                      %rows. Vertices are
                                                      %ordered
                                                      %counter-clockwise
                                                      %bottom then
                                                      %directly above's, eg
                                                      %a,b,c,a',b',c'.
                                                               
end
                                                


% Augment edges list by connecting vertices that
% are difference height 1 in up and down faces.
for i=1:n_sez-1
    e_slch = [e_slch ...
                     [(i-1)*n_esl+1:i*n_esl;...       %lower face nodes
                      i*n_esl+1:(i+1)*n_esl;...       %upper face nodes
                      e_sl(3,:)]...                 %store the arcs they
                     ];                             %may belong to.
                                                                       
end
n_pslch  = size(p_slch,2);                  %Update number of nodes, edges
n_eslch  = size(e_slch,2);                  %faces and prisms.
n_fslch  = size(f_slch,2);
n_prslch = size(pr_slch,2); 



%%%%%%%%%%%%%%%%%% Create n_chambers many meshed chambers %%%%%%%%%%%%%%%%%
% Interdiscal thickness is \nu*\epsilon_0
% Block disc thickness is \epsilon_0
% The first and last interdiscal spaces
% are supposed to have thickness 
% 1/2*\nu*\epsilon_0. This means that the
% z-axis H goes like .5I [C I] ... [C I] C .5I
% where the groupings are the (nu+1)*eps0 length
% If say there are n many C copies, there are 
% n-1 groupings then and n+1 interdiscal chambers
% where two are half-length. Ie 
% n = n_chambers - 1. 

% z-tic will give the base z-height of each
% I,C chamber in sequence. Careful though  
% that the first and last are only half as tall
% ie have height .5*nu*epsilon_0
disc = [1 0];
disc = repmat(disc,1,n_chambers-2);
disc = [0 disc 1 0];
disc = epsilon_0*disc;                  %Thickness of the C's

inter= [0 1];
inter= repmat(inter,1,n_chambers-2);
inter= [.5 inter 0 .5];
inter= nu*epsilon_0*inter;                     %Thickness of the I's.

z_tic= cumsum(inter+disc);                     %Compute right ends along z-
                                               %axis of the .5I {[C I]} C
                                               %.5I sequence.  
z_tic= z_tic - (inter+disc);                   %Compute left endpoint of 
                                               %intervals by subtracting
                                               %out their lengths.



p_3d  = p_ch;
e_3d  = e_ch;
f_3d  = f_ch;
pr_3d = pr_ch;

%vert_sl should be at p_3d(:,n_b4p-n_psl+1:n_psl))
%For adjoining slivers is p_3d(:,nb4p-npd+vert_sl)
%f_ch vertices in vert_sl should be swapped to 
%1st p_3d above.  All the rest should be adjusted
%by a n_b4p shift.  Because the call is given 
%relative end of vertices not beginning, the index
%rule conveniently won't need adjusting after each 
%gluein.

for i=1:2*(n_chambers-1)+1 
   disp(['Progress towards completing mesh: ' num2str(100*i/(...
        2*(n_chambers-1)+1 ...
        ))]);
   if     i == 1                        %First .5I chamber
             p_3d(1,:) = .5*nu*epsilon_0*p_3d(1,:); 
   elseif i == 2*(n_chambers-1)+1       %Final .5I chamber
       %Augment Vertex List
       newp = p_ch;                     %Don't duplicate the vertices are
       newp(1,:) = nu*epsilon_0*...     %glueing to at chamber's base. Note
                   newp(1,:);           %these are all vertices supplied by
       newp(:,vert_sl) = [];            %new chamber attachment, not just
                                        %at base. This moment determines
                                        %how the new vertex indices will
                                        %need to be shifted, and it depends
                                        %on their position. See below.

       n_b4p= size(p_3d,2);            %Number nodes in completed geometry                                 
       p_3d = [p_3d...                 %New nodes are translated and scaled
                z_tic(1,i)*[ones(1,n_pch-n_psl);zeros(4,n_pch-n_psl)]+...
                newp.*[.5*ones(1,n_pch-n_psl);ones(4,n_pch-n_psl)]...
               ];
                                       %.5 is only difference between this
                                       %and the general case.

       %Augment Faces List
       ram_f= f_ch;                    %Face structure for model chamber
       
       for j=1:n_psl                           %Renumber glued vertices
           overlap = find([f_ch(1:4,:);...
               zeros(1,size(f_ch,2))]==...
               vert_sl(j,1));
           ram_f(overlap) = n_b4p-n_psl+j;     %Set glued vertex value to
                                               %pre-existing index.
       end
       
                                                %Shift new index vertices
       shifts = zeros(1,n_pch);                 %Because sliver vertices
       shifts(vert_sl) = 1;                     %were scattered, the new 
       shifts = cumsum(shifts);                 %vertices aren't all 
                                                %index shifted by same
                                                %amount. -shifts(1,j) is
                                                %amount p_3d(:,j) is
                                                %shifted down in model
                                                %faces list.  Remains to
                                                %add n_b4p to get place in
                                                %p_3d
       overlap = ismember(f_ch(1:4,:),vert_sl); %Pre-existing nodes in face
       new     = 1-overlap;                     %New nodes in faces
       new     = new.*(f_ch(1:4,:)~=0);         %Only find the indices 
                                                %meaning vertices, not the
                                                %dummy 0 used to store a 
                                                %triangle face.
       new     = find([new;zeros(1,size(f_ch,2))]);
       for j=1:size(new,1)
           ram_f(new(j,1)) = f_ch(new(j,1)) -...%Make appropriate shift to
                             shifts(1,f_ch(new(j,1)))+... %match p_3d list.
                             n_b4p;
       end
                                                
       f_3d    = [f_3d ram_f];
       clear overlap new ram_f shifts
       
       %Augment prisms list
       ram_pr= pr_ch;                    %Face structure for model chamber
       
       for j=1:n_psl                           %Renumber glued vertices
           overlap = find([pr_ch(1:6,:);...
               zeros(1,size(pr_ch,2))]==...
               vert_sl(j,1));
           ram_pr(overlap) = n_b4p-n_psl+j;     %Set glued vertex value to
                                                %pre-existing index.
       end
       
                                                %Shift new index vertices
       shifts = zeros(1,n_pch);                 %Because sliver vertices
       shifts(vert_sl) = 1;                     %were scattered, the new 
       shifts = cumsum(shifts);                 %vertices aren't all 
                                                %index shifted by same
                                                %amount. -shifts(1,j) is
                                                %amount p_3d(:,j) is
                                                %shifted down in model
                                                %faces list.  Remains to
                                                %add n_b4p to get place in
                                                %p_3d
       overlap = ismember(pr_ch(1:6,:),vert_sl);%Pre-existing nodes in face
       new     = 1-overlap;                     %New nodes in faces
       new     = find([new;zeros(1,size(pr_ch,2))]);
       for j=1:size(new,1)
           ram_pr(new(j,1)) = pr_ch(new(j,1)) -...%Make appropriate shift to
                             shifts(1,pr_ch(new(j,1)))+... %match p_3d list.
                             n_b4p;
       end
                                                
       pr_3d    = [pr_3d ram_pr];
       clear overlap new ram_pr shifts       
   elseif mod(i,2) == 1 %Intermediate Chamber
                                                                            
       %Augment Vertex List
       newp = p_ch;                    %Don't duplicate the vertices are
       newp(1,:) = nu*epsilon_0*...     %glueing to at chamber's base. Note
                   newp(1,:);           %these are all vertices supplied by
       newp(:,vert_sl) = [];            %new chamber attachment, not just
                                        %at base. This moment determines
                                        %how the new vertex indices will
                                        %need to be shifted, and it depends
                                         %on their position. See below.
       
       n_b4p= size(p_3d,2);            %Number nodes in completed geometry                                 
       p_3d = [p_3d...                 %New nodes are translated and scaled
                z_tic(1,i)*[ones(1,n_pch-n_psl);zeros(4,n_pch-n_psl)]+...
                newp.*[ones(1,n_pch-n_psl);ones(4,n_pch-n_psl)]...
               ];

       %Augment Faces List
       ram_f= f_ch;                    %Face structure for model chamber
       
       for j=1:n_psl                           %Renumber glued vertices
           overlap = find([f_ch(1:4,:);...
               zeros(1,size(f_ch,2))]==...
               vert_sl(j,1));
           ram_f(overlap) = n_b4p-n_psl+j;     %Set glued vertex value to
                                               %pre-existing index.
       end
       
                                                %Shift new index vertices
       shifts = zeros(1,n_pch);                 %Because sliver vertices
       shifts(vert_sl) = 1;                     %were scattered, the new 
       shifts = cumsum(shifts);                 %vertices aren't all 
                                                %index shifted by same
                                                %amount. -shifts(1,j) is
                                                %amount p_3d(:,j) is
                                                %shifted down in model
                                                %faces list.  Remains to
                                                %add n_b4p to get place in
                                                %p_3d
       overlap = ismember(f_ch(1:4,:),vert_sl); %Pre-existing nodes in face
       new     = 1-overlap;                     %New nodes in faces
       new     = new.*(f_ch(1:4,:)~=0);         %Only find the indices 
                                                %meaning vertices, not the
                                                %dummy 0 used to store a 
                                                %triangle face.
       new     = find([new;zeros(1,size(f_ch,2))]);
       for j=1:size(new,1)
           ram_f(new(j,1)) = f_ch(new(j,1)) -...%Make appropriate shift to
                             shifts(1,f_ch(new(j,1)))+... %match p_3d list.
                             n_b4p;
       end
                                        
       f_3d    = [f_3d ram_f];
       clear overlap new ram_f shifts
       
       %Augment prisms list
       ram_pr= pr_ch;                    %Face structure for model chamber
       
       for j=1:n_psl                           %Renumber glued vertices
           overlap = find([pr_ch(1:6,:);...
               zeros(1,size(pr_ch,2))]==...
               vert_sl(j,1));
           ram_pr(overlap) = n_b4p-n_psl+j;    %Set glued vertex value to
                                               %pre-existing index.
       end
       
                                                %Shift new index vertices
       shifts = zeros(1,n_pch);                 %Because sliver vertices
       shifts(vert_sl) = 1;                     %were scattered, the new 
       shifts = cumsum(shifts);                 %vertices aren't all 
                                                %index shifted by same
                                                %amount. -shifts(1,j) is
                                                %amount p_3d(:,j) is
                                                %shifted down in model
                                                %faces list.  Remains to
                                                %add n_b4p to get place in
                                                %p_3d
       overlap = ismember(pr_ch(1:6,:),vert_sl);%Pre-existing nodes in face
       new     = 1-overlap;                     %New nodes in faces
       new     = find([new;zeros(1,size(pr_ch,2))]);
       for j=1:size(new,1)
           ram_pr(new(j,1)) = pr_ch(new(j,1)) -...%Make appropriate shift to
                             shifts(1,pr_ch(new(j,1)))+... %match p_3d list.
                             n_b4p;
       end
                                                
       pr_3d    = [pr_3d ram_pr];
       clear overlap new ram_pr shifts
   else                                 %Intermediate C Sliver Chamber                                                                                  
       %Augment Vertex List
       newp = p_slch;                   %Don't duplicate the vertices are
       newp(1,:) = epsilon_0*...        %glueing to at chamber's base. Note
                   newp(1,:);           %these are all vertices supplied by
       newp(:,1:n_psl) = [];            %new chamber attachment, not just
                                        %at base. This moment determines
                                        %how the new vertex indices will
                                        %need to be shifted, and it depends
                                        %on their position. See below.
       
       n_b4p= size(p_3d,2);            %Number nodes in completed geometry                                 
       p_3d =[p_3d...                  %New nodes are translated and scaled
              z_tic(1,i)*[ones(1,n_pslch-n_psl);zeros(4,n_pslch-n_psl)]+...
              newp.*[ones(1,n_pslch-n_psl);ones(4,n_pslch-n_psl)]...
             ];

       %Augment Faces List
       ram_f= f_slch;                  %Face structure for model sliver
       
       for j=1:n_psl                           %Renumber glued vertices
           overlap = find([f_slch(1:4,:);...   %In the sliver, glued
               zeros(1,size(f_slch,2))]==...   %vertices are first nodes
               j);
           ram_f(overlap) = n_b4p-n_pd+...
                             vert_sl(j,1);     %Set glued vertex value to
                                               %pre-existing index.
       end
       
                                                %Shift new index vertices
       new = find([f_slch(1:4,:);...            %Above n_psl iff we have 
                  zeros(1,n_fslch)]...          %it' s outside the sliver
                  >n_psl);                     
       ram_f(new) = f_slch(new)-n_psl+n_b4p;   %Because when joining a     
                                               %sliver all vert_sl are  
       f_3d    = [f_3d ram_f];                 %given first, the new
       clear overlap new ram_f                 %nodes are shifted the same
                                               %way

       %Augment prisms list
       ram_pr= pr_slch;                     %Prism structure for model sliver
       
       for j=1:n_psl                           %Renumber glued vertices
           overlap = find([pr_slch(1:6,:);...   %In the sliver, glued
               zeros(1,size(pr_slch,2))]==...   %vertices are first nodes
               j);
           ram_pr(overlap) = n_b4p-n_pd+...
                             vert_sl(j,1);     %Set glued vertex value to
                                               %pre-existing index.
       end
       
                                                %Shift new index vertices
       new = find([pr_slch(1:6,:);...            %Above n_psl iff we have 
                  zeros(1,n_prslch)]...          %it' s outside the sliver
                  >n_psl);                     
       ram_pr(new) = pr_slch(new)-n_psl+n_b4p;   %Because when joining a     
                                                 %sliver all vert_sl are  
       pr_3d    = [pr_3d ram_pr];                %given first, the new
       clear overlap new ram_pr                 %nodes are shifted the same
                                                 %way                                                                   
   end
end

%%%%%%%%%%%%% Adjust the nodes list so that mesh is for a %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%% true cone
p_3d(1:3,:) = [0 1 0;0 0 1;1 0 0]*p_3d(1:3,:); %Permute z to third position

for i=1:size(p_3d,2)
p_3d(1:3,i) = Fmap(p_3d(1:3,i));     %Map nodes to cone. Fmap is 
                                     %subroutine defined below.
end

save('Mesh','p_3d','f_3d','pr_3d');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
% plot mesh of prisms
% if plot_mesh

	% Mesh of prisms in the cone
    vert = p_3d(1:3,:);
    % Horizontal BOUNDARY faces of cone
    %figure
    %bd_chfaces = find(f_3d(5,:) == -1);
    %l=patch('Faces',f_3d(1:3,bd_chfaces)',...
    %    'Vertices',vert');
    %set(l,'facecolor','c')
	%axis equal
	%xlabel('x [\mu m]')
	%ylabel('y [\mu m]')
	%zlabel('z [\mu m]')
    %axis([-R_b*1.5 R_b*1.5 -R_b*1.5 R_b*1.5 0 H])
    % Exterior Lateral Chambers Cone
    %figure
    %bd_latfaces = find(...
    %    (f_3d(4,:)~=0).*(f_3d(5,:)==1)+...
    %     (f_3d(5,:)==-2));
    %l=patch('Faces',f_3d(1:4,bd_latfaces)',...
    %'Vertices',vert');
    %set(l,'facecolor','c')
	%axis equal
	%xlabel('x [\mu m]')
	%ylabel('y [\mu m]')
	%zlabel('z [\mu m]')
    %axis([-R_b*1.5 R_b*1.5 -R_b*1.5 R_b*1.5 0 H])
    % Outersliver boundary
    %figure
    %osliver = find(f_3d(5,:)==5);
    %l=patch('Faces',f_3d(1:4,osliver)',...
    %'Vertices',vert');
    %set(l,'facecolor','c')
	%axis equal
	%xlabel('x [\mu m]')
	%ylabel('y [\mu m]')
	%zlabel('z [\mu m]')
    %axis([-R_b*1.5 R_b*1.5 -R_b*1.5 R_b*1.5 0 H])
    % Sliver panels
    %figure
    %panel = find((f_3d(5,:) == 3)+...
    %             (f_3d(5,:) == 4));
    %l=patch('Faces',f_3d(1:4,panel)',...
    %'Vertices',vert');
    %set(l,'facecolor','c')
	%axis equal
	%xlabel('x [\mu m]')
	%ylabel('y [\mu m]')
	%zlabel('z [\mu m]')
    %axis([-R_b*1.5 R_b*1.5 -R_b*1.5 R_b*1.5 0 H])
    % Full cone plot
    %figure
    % bases                        
    %check = f_3d(4,:);                         %Extract row, says if have a 
    %bases = find(check == 0);                  %a face that is triangle or
    %l=patch('Faces',f_3d(1:3,bases)',...       %quadrilateral.
    %        'Vertices',vert');
    % lateral surface
    %quadrilateral = find(check ~= 0);
	%l=patch('Faces',f_3d(1:4,quadrilateral)',...
    %        'Vertices',vert'...
    %        );   
    %clear check
    %set(l,'facecolor','c')
	%axis equal
	%xlabel('x [\mu m]')
	%ylabel('y [\mu m]')
	%zlabel('z [\mu m]')
    %axis([-R_b*1.5 R_b*1.5 -R_b*1.5 R_b*1.5 0 H])
%%%%%
    %This below is to plot a single chamber
    %vert = [0 1 0;0 0 1;1 0 0]*p_slch(1:3,:);    %Permute z to third position 
	%l=patch('Faces',f_slch(1:3,1:n_sez*n_tsl)',...
%            'Vertices',vert');                 %With patch the polygons
                                               %need to be listed as rows
                                               %similarly each vertex needs
                                               %to be given as a row.
    %Use as a check to make sure spacing is as should be
    %max(abs(unique(vert(3,:))-...
    %nu*epsilon_0/(n_sez-1)*(cumsum(ones(1,n_sez))-1)))
    % lateral surface
	%l=patch('Faces',f_slch(1:4,n_sez*n_tsl+1:n_sez*n_tsl+(n_sez-1)*n_esl)',...
    %        'Vertices',vert'...
    %        );
	%set(l,'facecolor','c')
	%axis equal
	%xlabel('x [\mu m]')
	%ylabel('y [\mu m]')
	%zlabel('z [\mu m]')
    %axis([-R_b*1.1 R_b*1.1 -R_b*1.1 R_b*1.1 -.5 1.5])
	%view(3)
%end

end

function pmap = Fmap(p)
%Wrote code this way to make sure no 0/0 
%errors occurred when using logicals
%to piecewise define the map.

% scale nodes horizontally by 
% lambda(z) and the sliver nodes get shifted by 
% rule, which is x = Rte_{rho} + se_{rho} -->
% \lambda(z)Rte_{rho} + se_{rho} and equiv to
% x -- > x + [\lambda(z)-1]*Rte_{rho}

% Read geometric parameters
[ R_b,R_t,H,~,...
    ~] = data;
lambda = @(z) R_b/R_t - (R_b-R_t)/(H*R_t).*z;  %Function handle scaling
                                               %factor
pmap = zeros(3,1);
if p(1)^2+p(2)^2 <= R_t^2
    %In in cone, rescale cross sections
    pmap(1:2) = [lambda(p(3))*p(1);...
                 lambda(p(3))*p(2);];
    pmap(3)   = p(3);
else
    %In sliver do other map
    pmap(1:2) = [p(1) + (lambda(p(3))-1)*R_t*p(1)/sqrt(p(1)^2+p(2)^2);...
                 p(2) + (lambda(p(3))-1)*R_t*p(2)/sqrt(p(1)^2+p(2)^2)];
    pmap(3)   = p(3);
end
end

