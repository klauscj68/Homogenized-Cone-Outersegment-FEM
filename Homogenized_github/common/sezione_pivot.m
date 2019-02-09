% mesh del disco pivot, con e senza incisure
function [n_pd, n_int_pd, n_tri, n_sl, n_fo, ...
        sl2pd, fo2pd, pd2int_pd, p_pd, t_pd]=...
    sezione_pivot(R, taglia, n_ref, tol_R, theta_in, theta_fin, tol_angle)





% the local parameter n_ref:
% corresponds to n_ref_cyto when sezione_pivot is called by genemesh, for the generation of the mesh for the homogenized problem;
% corresponds to n_ref_id when sezione_pivot is called by genemesh, for the
% generation of the mesh for the diffusion of E* on the special disc

% counters
% n_pd is the number of nodes on the pivot disc
% n_int_pd is the number of nodes on the pivot disc not belonging to the
% sliver
% n_tri is the number of triangles of the pivot disc mesh
% n_sl is the number of nodes on the sliver

% correspondence matrices
% sl2pd contains the nodal indeces of nodes belonging to the sliver, in
% counterlock wise order.
% pd2int_pd associates to each node of the pivot disc: 
%   zero, if the node is on the sliver; 
%   a progressive number counting nodes belonging to the interior of the disc.

% mesh
% p_pd nodal coordinates vector
% t_pd adjacence matrix of the triagles



% genera la mesh non incisa
fun='disc_sliver';
[p_pd,e_pd,t_pd]=initmesh(fun,'Hmax',taglia); 
% n_ref raffittimenti regolari
for j=1:n_ref
    [p_pd,e_pd,t_pd]=refinemesh(fun,p_pd,e_pd,t_pd,'regular'); 
end

% number of nodes on the pivot disc
n_pd=size(p_pd,2);
% number of triangles on the pivot disc
n_tri=size(t_pd,2);

% nodes on the sliver and the folds
% choose all nodes at a distance >=R-tol_R from the origin
R2pd=find(abs(p_pd(1,:)+1i*p_pd(2,:))>=R-tol_R);



% compute the corresponding angles
angoli=angle(p_pd(1,R2pd)+1i*p_pd(2,R2pd));
% trasporta da (-pi,pi] a [0,2*pi)
negativi=find(angoli<0-tol_angle);
angoli(negativi)=angoli(negativi)+2*pi;

% select the nodes belonging to the sliver 
ang_ind=angoli>=theta_in-tol_angle & angoli<=theta_fin+tol_angle;
sl2pd=R2pd(ang_ind);
% number of nodes on the sliver
n_sl=length(sl2pd);
% ordina nel senso delle anomalie crescenti
[~,I]=sort(angoli(ang_ind));
sl2pd=sl2pd(I);

% select the nodes belonging to the folds
angoli(angoli<theta_in+tol_angle)=angoli(angoli<theta_in+tol_angle)+2*pi;
ang_ind=angoli>=theta_fin-tol_angle;
fo2pd=R2pd(ang_ind);
% number of nodes on the sliver
n_fo=length(fo2pd);
% ordina nel senso delle anomalie crescenti
[~,I]=sort(angoli(ang_ind));
fo2pd=fo2pd(I);




% determines pd2int_pd

% the negation of ismember has 0 in correspondence of the nodes on the
% sliver, and 1 in correspondence of the remaining nodes
interni=~ismember(1:n_pd,sl2pd);
% number of the nodes which are not on the sliver (equal to the number of ones in the vector interni)
n_int_pd=sum(interni);
pd2int_pd=cumsum(interni).*interni;



% get read of the fourth line which is the zone index
t_pd=t_pd(1:3,:);


return
