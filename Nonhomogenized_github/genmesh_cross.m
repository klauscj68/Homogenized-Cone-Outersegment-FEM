function [ p_pd,e_pd,t_pd,e_arcpd,n_pd,n_ed,n_tri,n_arced,...
           tri_disc,tri_sl,...
           edgarc1,edgarc2,edgarc3,edgarc4,edgarc5,...
           vert_disc,vert_sl,...
           vertarc1,vertarc2,vertarc3,vertarc4,vertarc5]...
           = genmesh_cross(R_b,R_t,H,theta_in,theta_fin,...
                           epsilon_0,nu,sigma,...
                           taglia,tol_R,tol_angle,n_ref)
%Subroutine of genmesh that handles meshing the cross section.
%   Mesh a cone cross section before scaling and output all necessary lists
%   used to build the 3d geometry. This is a counterpart to sezione_pivot.
%   The outputs p_pd has been augmented to include an additional 2 rows
%   that record membership in the volume disc, volume sliver in the first
%   row by recording a 1 or a 2, and in the second row by recording a 1 
%   through 5.  Note that nodes belonging to interfaces are recorded as 
%   belonging to the sliver or the index ordered last arc it belongs.  (In
%   other words prior membership is overwritten.)

%   e_pd's first two rows are the vertex index number recorded from least
%   to greatest. Its third row records which domain arc it may belong
%   (indexed 1-5) or no membership (index 0).

%   e_arcpd are those edges belonging specifically to the
%   domain arcs.  Similarly first two rows are its vertex indices from
%   least to greatest.  Its third row records its edge index in e_pd. Its
%   4th row records which arc (1 to 5) that it belongs.  

%p_pd are nodal points but e_pd are only the edges on the boundary arcs not
%all edges in the mesh. t_pd are otherwise all mesh triangles.
fun='pdegeom_disc_sliver';
[p_pd,e_arcpd,t_pd]=initmesh(fun,'Hmax',taglia); %Hmin gives minimum length

% Refine initial disc mesh specified number of times
for j=1:n_ref
    [p_pd,e_arcpd,t_pd]=refinemesh(fun,p_pd,e_arcpd,t_pd,'regular'); 
end

% number of nodes on disc cross section
n_pd=size(p_pd,2);
% number of arc edges on disc cross section
n_arced=size(e_arcpd,2);
% number of triangles on disc cross section
n_tri=size(t_pd,2);

% Build the actual edges list and the number
% of edges.
e_pd = zeros(2,3*n_tri);       %At most the triangles are disjoint
n_ed = 0;                      %Initialize edge count

for i=1:n_tri
    triangle = sort(t_pd(1:3,i));       %Order vertices in ith triangle
                                        %from least to greatest.
                                        
    triangle = [triangle(1,1) triangle(2,1);...         %Store the edges 
                triangle(2,1) triangle(3,1);...         %comprising the
                triangle(1,1) triangle(3,1)]';          %triangle by its 
                                                        %edges which are
                                                        %stored by their
                                                        %vertices that are 
                                                        %least to greatest
                                                        %themselves. Each
                                                        %edge is stored by
                                                        %least to greatest
                                                        %vertex index.
    if i==1
        e_pd(:,1:3) = triangle;                         %First triangle 
                                                        %will not duplicate
                                                        
        n_ed        = 3;                                %Now are 3 edges
    else
        for j=1:3                                       %Cycle over the 
                                                        %edges of the given
                                                        %triangle
            
            current     = 1;                            %Initialize counter
                                                        %for current
                                                        %position in edges
                                                        %list.
                                                        
            found       = 0;                            %Initialize logical
                                                        %to see if current
                                                        %edge has been
                                                        %found.
        while (found == 0) && (current <= n_ed)
            if (triangle(1,j)==e_pd(1,current)) &&...
                    (triangle(2,j)==e_pd(2,current))
                found = 1;                              %Record if found.
            else
                current = current + 1;                  %If not, cycle to 
                                                        %the next.
            end
        end
        if found == 0
            e_pd(:,current) = triangle(:,j);            %If never found 
                                                        %record it in the
                                                        %edge list as a new
                                                        %edge
                                                        
            n_ed = n_ed + 1;                            %Update total 
                                                        %number of edges
                                                        %count.
        end
        end
    end
end
if n_ed < 3*n_tri
   e_pd = e_pd(:,1:n_ed);                               %If necessary 
                                                        %delete the padded
                                                        %zeros.
end


% Write into the e_arcpd list what edge number
% these boundary arc edges correspond to.
e_arcpd(1:2,:) = sort(e_arcpd(1:2,:));              %Write the edges in 
                                                    %their normal form.
                                                    
e_arcpd(3,:)   = zeros(1,n_arced);                  %Overwrite the 
e_arcpd(6:7,:) = [];                                %arclength 
                                                    %parameterization
                                                    %terms and also the
                                                    %domains lying to the
                                                    %left and right of
                                                    %edge information.
e_arcpd(4,:)   = [];

for i=1:n_arced
    %Cycle over all boundary arcs until they
    %are found in e_pd
    found = 0;
    count = 1;
    while (found == 0) && (count <=n_ed)
        if (e_arcpd(1,i)==e_pd(1,count)) &&...     %If you found it, record
                (e_arcpd(2,i)==e_pd(2,count))      %it
            found = 1;
            e_arcpd(3,i) = count;                  %Record what's actual
                                                   %edge number.
        else
            count = count + 1;
        end
    end
end

% Add in a third row to e_pd to record which arc
% if any this edge belongs.
e_pd = [e_pd; zeros(1,n_ed)];
for i=1:5
    check = e_arcpd(4,:);               %4th row records arc membership
    found = find(check == i);           %Find columns belong to ith arc
    index = e_arcpd(3,found);           %Extract edge indices on arc
    e_pd(3,index) = i;                  %These edges record arc membership
end
    
%%%%% Find the triangles belonging to the sliver
% After fact realized the quicker routine to find
% is reduce it to single row vector and then use 
% 'find' because that returns the column then 
% which is ultimately what after.  I keep it this
% way only to demonstrate that ind2sub command
% exists
tri_sl = t_pd.*([zeros(3,n_tri);ones(1,n_tri)]);   %Have nonzeros
                                                   %only in last row so
                                                   %that can retain (m,n)
                                                   %position in 'IND2SUB'
                                                   %call but so that call
                                                   %of 'find' only grabs
                                                   %domains not vertex
                                                   %indices.

tri_sl = find(tri_sl == 2);                        %Find triangles in the 
                                                   %sliver part of domain.
                                                   
[~,tri_sl] = ind2sub([4 n_tri],tri_sl);            %Extract those columns 
                                                   %associated to the
                                                   %triangle. Output is as
                                                   %column vector.
                                                   

%%%%% Find the triangles belonging to the disc
tri_disc = t_pd.*([zeros(3,n_tri);ones(1,n_tri)]); %Have nonzeros
                                                   %only in last row so
                                                   %that can retain (m,n)
                                                   %position in 'IND2SUB'
                                                   %call but so that call
                                                   %of 'find' only grabs
                                                   %domains not vertex
                                                   %indices.

tri_disc = find(tri_disc == 1);                    %Find triangles in the 
                                                   %sliver part of domain.
                                                   
[~,tri_disc] = ind2sub([4 n_tri],tri_disc);        %Extract those columns 
                                                   %associated to the
                                                   %triangle. Output is as
                                                   %column vector.
%%%%% Find the edges belonging to the geometric 
%%%%% arcs of the domain

edgarc1 =  e_arcpd(3,find(e_arcpd(4,:)==1));
edgarc2 =  e_arcpd(3,find(e_arcpd(4,:)==2));
edgarc3 =  e_arcpd(3,find(e_arcpd(4,:)==3));
edgarc4 =  e_arcpd(3,find(e_arcpd(4,:)==4));
edgarc5 =  e_arcpd(3,find(e_arcpd(4,:)==5));
                                                   
                                                   
%%%%%% Find nodes in \bar{sliver} volume by finding 
%%%%%% vertices associated to the found sliver 
%%%%%% triangles.

vert_sl = t_pd(1:3,tri_sl);                        %Extract sliver 
                                                   %triangles indexed by
                                                   %their vertices. Each
                                                   %entry is the index of a
                                                   %vertex belonging to
                                                   %sliver. Only issue is
                                                   %that there may be
                                                   %duplicates.

vert_sl = reshape(vert_sl,3*size(tri_sl,1),1);     %Reshape it into a 
                                                   %vector.
                                                   
vert_sl = unique(vert_sl);                         %Parse down any 
                                                   %duplicates
                                                   
%%%%%% Find nodes in \bar{disc} volume by finding 
%%%%%% vertices associated to the found disc 
%%%%%% triangles.

vert_disc = t_pd(1:3,tri_disc);                    %Extract sliver 
                                                   %triangles indexed by
                                                   %their vertices. Each
                                                   %entry is the index of a
                                                   %vertex belonging to
                                                   %sliver. Only issue is
                                                   %that there may be
                                                   %duplicates.
                                                   
vert_disc = reshape(vert_disc,3*size(tri_disc,1),1);%Reshape it into a 
                                                    %vector.
                                                   
vert_disc = unique(vert_disc);                      %Parse down any 
                                                    %duplicates
                                                                  
%%%%% Find the nodes belonging to the arcs

%Nodes of Arc 1
vertarc1 = e_pd(1:2,edgarc1);                      %Find columns of 
                                                   %e_arcpd with this arc
                                                   %and then extract the 
                                                   %edge column index of
                                                   %these. To extract their
                                                   %vertices.  Probably
                                                   %could make faster by
                                                   %going direct through
                                                   %first two rows of
                                                   %e_arcpd instead.
       
vertarc1 = reshape(vertarc1,2*size(edgarc1,2),1);
vertarc1 = unique(vertarc1);                       %Convert vertices 
                                                   %belonging to these
                                                   %vectors into column
                                                   %vector and eliminate
                                                   %any duplicates.

%Nodes of Arc 2
vertarc2 = e_pd(1:2,edgarc2);                      %Find vertex at these 
                                                   %edges.
vertarc2 = reshape(vertarc2,2*size(edgarc2,2),1);
vertarc2 = unique(vertarc2);                       %Convert vertices 
                                                   %belonging to these
                                                   %vectors into column
                                                   %vector and eliminate
                                                   %any duplicates.

%Nodes of Arc 3
vertarc3 = e_pd(1:2,edgarc3);                      %Find vertex at these 
                                                   %edges.
vertarc3 = reshape(vertarc3,2*size(edgarc3,2),1);
vertarc3 = unique(vertarc3);                       %Convert vertices 
                                                   %belonging to these
                                                   %vectors into column
                                                   %vector and eliminate
                                                   %any duplicates.
                                                   
%Nodes of Arc 4
vertarc4 = e_pd(1:2,edgarc4);                      %Find vertex at these 
                                                   %edges.
vertarc4 = reshape(vertarc4,2*size(edgarc4,2),1);
vertarc4 = unique(vertarc4);                       %Convert vertices 
                                                   %belonging to these
                                                   %vectors into column
                                                   %vector and eliminate
                                                   %any duplicates.
                                                   
%Nodes of Arc 5
vertarc5 = e_pd(1:2,edgarc5);                      %Find vertex at these 
                                                   %edges.
vertarc5 = reshape(vertarc5,2*size(edgarc5,2),1);
vertarc5 = unique(vertarc5);                       %Convert vertices 
                                                   %belonging to these
                                                   %vectors into column
                                                   %vector and eliminate
                                                   %any duplicates.
                                                   
%%%%% Augment the original p_pd, e_pd, lists
%%%%% so to record membership in volume disc, 
%%%%% volume sliver, arc1 through 5.  This is done
%%%%% by adding 7 new rows with 0 and 1 used to 
%%%%% indicate membership. Standardize the t_pd
%%%%% list to have last two rods record this
%%%%% membership.

% Update p_pd
% Note overlaps will only belong to one of these
ram = zeros(2,n_pd);
ram(1,vert_disc') = 1; %Record 1 in nodes in volume disc
ram(1,vert_sl')   = 2; %Record 2 in nodes in volume sliver
ram(2,vertarc1')  = 1; %Record nodes in arc 1
ram(2,vertarc2')  = 2; %Record nodes in arc 2
ram(2,vertarc3')  = 3; %Record nodes in arc 3
ram(2,vertarc4')  = 4; %Record nodes in arc 4
ram(2,vertarc5')  = 5; %Record nodes in arc 5
p_pd = [p_pd;ram];     %Write into p_pd

end

