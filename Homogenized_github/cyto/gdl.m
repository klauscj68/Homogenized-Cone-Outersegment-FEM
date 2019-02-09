function [n_gdl, gdl_vol, gdl_sd, gdl_sl, gdl_fo]=...
    gdl(n_pd, n_int_pd, n_sl, n_fo, ...
    n_sez, n_sd, ...
    sl2pd, fo2pd, pd2int_pd, sd2sez, flag_model, inspect)

% degrees of freedom (dof's) corresponding to the volume, the special discs and the
% sliver

% order of the dof's in the global matrices
% 3D model (flag_model=1)
%   first nodes in the volume, equal to the number of sections (n_sez) times the nodes on the pivot disc (n_pd)
%   then the nodes of the special discs not belonging to the sliver (n_int_pd*n_sd), 
%   not additional dof must be considered for the sliver, since the
%   solution on the sliver is obtained as a trace of the solution in the
%   volume
%  transversal well-stirred model (flag_model=2)
%   1 gdl for all the nodes belonging to a cross section;
%   the dof's are ordered on the sense of the increasing z
%  global well-stirred model (flag_model=3)
%   one dof for all the nodes

% counters
% n_pd is the number of nodes on the pivot disc
% n_int_pd is the number of nodes of the pivot disc not belonging to the sliver
% n_tri is the number of triangles on the pivot disc
% n_sl is the number of nodes on the sliver
% n_sez is the number of tranversal sections in the volume
% n_sd is the number of special discs

% correspondences
% sl2pd contains the indeces of the nodes of the trace of the sliver in the
% pivot disc, in cunterclock wise, in the ordering of the pivot disc
% pd2int_pd associates at each node on the pivot disc: 
%   zero, if the node is on the sliver; 
%   a progressive number counting othe nodes of the pivot disc
%   not on the sliver trace
% sd2sez associates to each special disc the section number superimposed


% output
% n_gdl is the total number of dof's for each variable (u e v)
% the total number of dof's of the problem is equal to 2*n_gdl
% gdl_vol(j,k)    is the dof of the node of index j in the numbering of the pivot disc, in section k
% gdl_sd(j,k)     is the dof of the node of index j in the numbering of the pivot disc, in the special disc k
% gdl_os(j,k)     is the dof of the node of index j in the sliver (in the numbering of the sliver), in section k

fprintf('\nDetermine correspondence between nodes and degrees of freedom\n');


switch flag_model
    case 1
        % modello 3D

        % dof couter
        n_gdl=0;

        % dof's of nodes in the volume
        gdl_vol=n_pd*repmat(0:n_sez-1,n_pd,1)+ repmat((1:n_pd)',1,n_sez);
        % update dof's
        n_gdl=n_gdl+n_sez*n_pd;

        % dof's of nodes in the special discs
        gdl_sd=zeros(n_pd,n_sd);
        for k=1:n_sd
            % for nodes on the special discs, where pd2int_pd is zero, take the
            % dof of the corresponding node in the volume
            % for each internal node, where pd2int_pd is not zero, creates a new dof which reads from pd2int_pd
            % dof if all nodes of the special disc coincide with the corresponding nodes in the volume 
            gdl_c=n_pd*(sd2sez(k)-1)+ (1:n_pd);
            % dof's belonging to internal nodes which are different from
            % nodes of the volume
            gdl_nc=n_gdl+pd2int_pd;
            % combine in suitable way:
            % in positions where pd2int_pd is zero (node on the sliver common to corresponding node in the volume), takes the dof from gdl_c;
            % in positions where pd2int_pd is non zero (new dof), take the dof from gdl_nc;
            gdl_sd(:,k)=(gdl_c.*(pd2int_pd==0)+gdl_nc.*(pd2int_pd~=0))';
            % aggiorna n_gdl, sommando quelli appena aggiunti
            n_gdl=n_gdl+n_int_pd;
        end

        % dof in the sliver, no new dof is created
        gdl_sl=n_pd*repmat(0:n_sez-1,n_sl,1)+repmat(sl2pd',1,n_sez);

        % dof in the folds, no new dof is created
        gdl_fo=n_pd*repmat(0:n_sez-1,n_fo,1)+repmat(fo2pd',1,n_sez);
        
        
    case 2
        % modello well-stirred trasversale

        % i gdl sono pari al numero delle sezioni
        n_gdl=n_sez;
    
        % i nodi di una sezione corrispondono allo stesso gdl
        gdl_vol=repmat(1:n_sez,n_pd,1);

        % gdl dei nodi dei dischi speciali
        % i nodi di ciascun disco speciale hanno il gdl
        % corrispondente ai nodi sulla circonferenza in comune con l'outer shell
        % che a sua volta è quello dei nodi dell'interior alla stessa quota
        gdl_sd=repmat(sd2sez,n_pd,1);

        % gdl dei nodi nell'outer shell
        gdl_sl=repmat(1:n_sez,n_sl,1);



    case 3
        % modello well-stirred globale

        % un solo gdl per tutti i nodi
        n_gdl=1;

        % volume
        gdl_vol=ones(n_pd,n_sez);
        
        % dischi speciali
        gdl_sd=ones(n_pd,n_sd);
        
        % outer shell
        gdl_sl=ones(n_sl,n_sez);

        
end % switch

% % stampe di controllo
% n_gdl, gdl_vol, gdl_sd, gdl_os
% for k=1:n_inc
%     'incisura',k
%     gdl_inc{k}
% end

% informativa
fprintf('\nGeneral mesh: %i overall dof''s (u e v)\n', 2*n_gdl);
presskey(inspect);

return
