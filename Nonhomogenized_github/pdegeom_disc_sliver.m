function [ x,y ] = pdegeom_disc_sliver( bs,s )
%Geometry file used by initmesh to encode meshed cross section w/sliver
%   In the help file of initmesh, it offers that geometries may be encoded
%   using a function in the template of 'pdegeom.' For example call help
%   pdegeom to see this template.  Note that bs and s are internal values
%   set by matlab's call of initmesh and not specified by the user.  Also
%   note that 'pdegeom' is set to call a variable number of inputs.  The
%   'nargin' command lets you recognize that number. The geometric cross
%   section is a disc (region 1) with an added radial sliver (region 2)
%   immersed in R2 (region 0). Regions that aren't to be meshed should be
%   given the number 0 according to matlab documentation. Rmbr want to 
%   build global mesh so that  light comes from the top which is the 
%   smaller end.

%   For more information, see this website
%   https://www.mathworks.com/help/pde/ug/create-geometry-using-a-geometry-function.html?requestedDomain=www.mathworks.com

%   Code to plot this geomtry file but should be called
%   outside of this function.
%   1. Gives geometry without mesh
%   pdegplot(@pdegeom_disc_sliver,'EdgeLabels','on','SubdomainLabels','on')
%   axis equal
%   2. Generates lists.  Note that the triangles 't' list has 4 rows
%      whose first three rows give vertex indices (oriented) and whose
%      final row indicates the region of pdegeom (ie 1 or 2) that the
%      triangle belongs.
%   [p,e,t] = initmesh(@pdegeom_disc_sliver); ",'Hmax',Number) restricts mesh size
%   3. Plot the geometry with mesh
%   pdeplot(p,e,t)
%   title('Mesh Plot')
%   axis equal

%Make a call to the geometric data defining the problem.
[~,R_t,~,theta_in,theta_fin,epsilon_0,~,sigma,~,...
    ~] = data();
%R_t       =  1;
%sigma     =  1;
%epsilon_0  = .1;

%According to matlab the output depends on the number of arguments 
%entered.
switch nargin
    case 0
        nbs = 5;  %Number of boundary segments is two circular 
                  %arcs for the inner disc, two straight edge
                  %radial extensions, and one final circular 
                  %arc. Enclosing the outer sliver.
        x = nbs;
        return
    case 1
        %One input is bs. 
        %"Here bs is a vector of edge segments. Your 
        %function returns d as a matrix with one column for each edge 
        %segment specified in bs. This is different than in case 2 FYI. 
        %The rows of D are:
        %Start parameter value
        %End parameter value
        %Left region label, where "left" is with respect to the direction from the start to end parameter value
        %Right region label. Region label is the same as subdomain number. 
        %The region label of the exterior of the geometry is 0." 
        nbs = 5;                   %Number of boundary edge segments
        
        D = zeros(4,nbs);          %For storing this mesh edge 
                                   %'adjacency' matrix as 
                                   %described in matlab. Each
                                   %column is associated to an
                                   %edge.  The first two rows
                                   %are abstract parameters (that
                                   %give the range of params) 
                                   %The first two rows are the start 
                                   %and enpoint of an interval 
                                   %parameterization for the intended 
                                   %boundary edge. According to matlab
                                   %documentation. The boundary is best 
                                   %if param'd by an affine function 
                                   %of its arclength. That specific param
                                   %isn't given until the nargin == 2
                                   %case where how they map to
                                   %boundary points is defined.
                                   %That is the parameterization 
                                   %associated to these isn't given
                                   %until the final case call of 
                                   %pdegeom. The last two rows indicate
                                   %the region to left and right
                                   %of the edge in the from above
                                   %counterclockwise sense.
                  
        %D(:,1) = [0 pi 1 0]';      %First boundary edge is the 
                                   %upper part of disc parameterized 
                                   %by the angle theta
        D(:,1) = [theta_fin 2*pi+theta_in 1 0]';                          

        %D(:,2) = [-pi 0 1 2]';     %Second boundary edge is the
                                   %lower part of disc parameterized
                                   %by the angle theta
        D(:,2) = [theta_in theta_fin 1 2]';
                       
        D(:,3) = [0 1 2 0]';       %Third boundary edge is the
                                   %left-side radial straight edge
                                   %which has thickness the sliver
                                   %and that is parameterized by the 
                                   %unit interval
                      
        D(:,4) = [0 1 2 0]';       %Fourth boundary edge is the 
                                   %right-side radial straight edge
                                   %which has thickness the sliver
                                   %and that is parameterized by the
                                   %unit interval.
                     
        %D(:,5) = [-pi 0 2 0]';    %Fifth boundary edge is the bottom
                                   %most sliver lip that has radius 
                                   %r + sliver thickness
        D(:,5) = [theta_in theta_fin 2 0]';
        
        x = D(:,bs);               %bs is a vector of indices matlab
                                   %will call into this function.  
                                   %D(:,bs) now gives all columns of D
                                   %which match the indices entered by
                                   %bs.
        return
    case 2
        %Here the inputs called by matlab are bs and s.
        %"s is an array of arc lengths, and bs is a scalar or an array the 
        %same size as s giving edge numbers. If bs is a scalar, then it 
        %applies to every element in s. Your function returns x and y, 
        %which are the x and y coordinates of the edge segments specified 
        %in bs at parameter value s. The x and y arrays have the same size 
        %as s."
        %Believe the parameterizations should be written to coincide with
        %above for which regions are on the left and which are on the right
        x     = zeros(size(s));
        y     = zeros(size(s));
        
        if numel(bs) == 1          %bs might need scalar expansion
            bs = bs*ones(size(s)); %expand bs to size(s) array. 
        end                        %Presumambly this case has
                                   %all the s values corresponding
                                   %to same edge in matlab's call.
        
        %The use of "find( )" will allow us to proceed edge by edge
        %Note that matlab can call entries of an array using only 
        %a single index by cycling down successive columns. In
        %particular "find( )" outputs these 1d indices and we'll
        %use those to store the appropriate values of x and y for
        %the intended boundary point.
        
        %First edge and second edge on same circle
        found        = find(bs < 3);
        x(found)     = R_t*cos(s(found));
        y(found)     = R_t*sin(s(found));
        
        %Third edge is the sliver line on left
        found       = find(bs == 3);
        %x(found)    = -R_t - sigma*epsilon_0*s(found);
        %y(found)    = 0;
        x(found)    = (R_t+sigma*epsilon_0*s(found))*cos(theta_in);
        y(found)    = (R_t+sigma*epsilon_0*s(found))*sin(theta_in);
        
        %Fourth edge is the sliver line on right
        found       = find(bs == 4);
        %x(found)    = R_t + sigma*epsilon_0*(1-s(found));           %Should this be 1-s for orientation not just s(found)?
        %y(found)    = 0;
        x(found)    = (R_t+sigma*epsilon_0*(1-s(found)))*cos(theta_fin);
        y(found)    = (R_t+sigma*epsilon_0*(1-s(found)))*sin(theta_fin);
        
        
        %Fifth edge is bottom sliver circle part
        found       = find(bs == 5);
        x(found)    = (R_t + sigma*epsilon_0)*cos(s(found));
        y(found)    = (R_t + sigma*epsilon_0)*sin(s(found));

        return
end
end

