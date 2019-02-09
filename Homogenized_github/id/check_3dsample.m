function [ V_cone,V_sl,p_sample_cone,p_sample_sl,...
           bbest_cone, bbest_sl ] =...
           check_3dsample(p_3d,pr,tri,p_cross,Sp,...
                          rinc,thetainc,nzsec,chinc,...
                          n_sd,Z_sd,sddepth)
%Master script to eval spline at all cone/sl samples given by incs
%   TO save time, pr should only have the 6 rows for vertex indices.
%   tri=pr(1:3,:) should be prebuilt. p_cross = p(1:2,:) should be
%   prebuilt.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the data
R_b = 3;
R_t = 1;
H   = 15;
nu  = 1;
n_chambers = 100;
epsilon_0 = H/((1+nu)*(n_chambers-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%p_3d was given in the physical cone space, but check_find pr needs points
%to be in the uniform R_t cylinder setting
p_3d = p_3d(1:3,:);
p_3d = check_3drevertp( p_3d );

%Getting ready to check if need to rotate cylinder points into
%common [-pi/2,pi/2]
geo_type = input('Is this Nonhomogenized(0) or Homogenized(1):');
R = [0 1 0;... cos(theta) -sin(theta) 0
         -1 0 0;...sin(theta) cos(theta)  0
         0  0 1];%     0          0       1
if geo_type == 1
    %   Rotate homogenized geometry to be coincident with Nonhomogenized
    %   In nonhom, sliver was [-pi/2,pi/2] while in hom was [0,pi] so we
    %   rotate by angle -pi/2
    p_3d = R*p_3d;
end

%Generate the points where we take samples
[p_sample_cone,p_sample_sl,~] =...
    check_3dgenpts(rinc,thetainc,nzsec,chinc,...
                                       n_sd,Z_sd,sddepth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluate the spline at all sample points
%if geo_type == 1 %Homogenized
%    p_sample = [p_sample_cone;
%                zeros(1,size(p_sample_cone,2))]; %pad extra zeros because p_sample_sl now in 4th gives its hom gdl index
%    ram = [R*p_sample_sl{1}(1:3,:);p_sample_sl{1}(4,:)]; %Rotate the hom's sliver nodes back to nonhom position
%    p_sample = [p_sample ram];
%else             %Nonhomogenized
%    p_sample = [p_sample_cone;...
%                zeros(1,size(p_sample_cone,2))]; %pad extra zeros because p_sample_sl now in 4th gives its hom gdl index
%   for i=1:size(p_sample_sl,2)
%       ram = [R*p_sample_sl{i}(1:3,:);p_sample_sl{i}(4,:)]; %Rotate the hom's sliver nodes back to nonhom position
%       p_sample = [p_sample ram];
%   end
%end

if geo_type == 1 %Handle the homogenized case
    V_cone = zeros(1,size(p_sample_cone,2));
    V_sl   = zeros(1,size(p_sample_sl{1},2));
    
    %Compute values at the internal cone nodes
    bbest_cone = zeros(1,size(V_cone,2));
    [~,~,~,pt_min,pt_max] = check_3dfindpr(p_3d,pr,tri,p_cross,...
                                           p_sample_cone(:,1),.01);
    wbar = waitbar(0,'Progress Through Volume Samples');
    for i=1:size(V_cone,2)
        %Find prism containing point
        [prfound,bbest_cone(i),~,~,~]=...
            check_3dfindpr(p_3d,pr,tri,p_cross,...
                           p_sample_cone(:,i),.01,pt_min,pt_max);
        %Extract nodes of that prism and value of spline at those points
        nodes = pr(1:6,prfound);
        nodals= Sp(nodes);
        %Evaluate the spline at the sample point
        P = p_3d(:,nodes);
        %Check that pr is formatted as expected
        test1 = unique(P(3,1:3));
        test2 = unique(P(3,4:6));
        if (size(test1,2)~=1)||(size(test2,2)~=1)
            error([num2str(prfound) '^th prism is not ordered with its base ' ... 
                   'the first 3 entries and cap face the last 3 entries '...
                   'pr list.']);
        end
        %   Have the system
        %   x = P[(1-eta)I3 0;0 eta*I3][I3;I3]\varphi
        %   s = nodals*[(1-eta)I3 0;0 eta*I3][I3;I3]\varphi
        %   varphi = [1;0;0]+[-1 -1;1 0;0 1]*[xi1;xi2]
        %   Goal: From x compute s
        eta = (p_sample_cone(3,i)-P(3,1))/...
              (P(3,4)-P(3,1)); %z-value interpolated between extremes
        M = P*...
            [(1-eta)*eye(3) zeros(3,3);...
            zeros(3,3)  eta*eye(3)]*...
            [eye(3);eye(3)];
        %Holds x = M*varphi = Me_{1} + M*[-1 -1;1 0;0 1]*[xi1;xi2]
        x = p_sample_cone(:,i);
        x = x - M*[1;0;0];
        M = M*[-1 -1;1 0;0 1];
        M_LI = (M'*M)^(-1)*M';
        xi   = M_LI*x; clear x;
        Lag  = [xi;eta];
        V_cone(i) = nodals*...
                [(1-Lag(3))*eye(3) zeros(3,3);...
                   zeros(3,3) Lag(3)*eye(3)]*...
                [eye(3);eye(3)]*...
                [1-Lag(1)-Lag(2);Lag(1);Lag(2)];
    
        %Update the waitbar
        waitbar(i/size(V_cone,2),wbar)
    end
    close(wbar)
    
    %Rotate sliver samples back to their [-pi/2,pi/2] position
    %but preserve the 4th row gdl index
    p_sample_sl{1} = [R zeros(3,1); zeros(1,3) 1]*p_sample_sl{1};
    for i=1:size(V_sl,2)
        %4th row of p_sample_sl stored its gdl index
        V_sl(i) = Sp(p_sample_sl{1}(4,i));
    end
    bbest_sl = NaN;
else
    %Nonhomogenous Case 
    V_cone = zeros(1,size(p_sample_cone,2));
    V_sl   = cell(1,size(p_sample_sl,2));
    for i=1:size(p_sample_sl,2)
        V_sl{i} = zeros(1,size(p_sample_sl{i},2));
    end
    
    %Compute values at the internal cone nodes
    bbest_cone = zeros(1,size(V_cone,2));
    [~,~,~,pt_min,pt_max] = check_3dfindpr(p_3d,pr,tri,p_cross,...
                                           p_sample_cone(:,1),.01);
    wbar = waitbar(0,'Progress Through Volume Samples');
    for i=1:size(V_cone,2)
        %Find prism containing point
        [prfound,bbest_cone(i),~,~,~]=...
            check_3dfindpr(p_3d,pr,tri,p_cross,...
                           p_sample_cone(:,i),.01,pt_min,pt_max);
        %Extract nodes of that prism and value of spline at those points
        nodes = pr(1:6,prfound);
        nodals= Sp(nodes);
        %Evaluate the spline at the sample point
        P = p_3d(:,nodes);
        %Check that pr is formatted as expected
        test1 = unique(P(3,1:3));
        test2 = unique(P(3,4:6));
        if (size(test1,2)~=1)||(size(test2,2)~=1)
            error([num2str(prfound) '^th prism is not ordered with its base ' ... 
                   'the first 3 entries and cap face the last 3 entries '...
                   'pr list.']);
        end
        %   Have the system
        %   x = P[(1-eta)I3 0;0 eta*I3][I3;I3]\varphi
        %   s = nodals*[(1-eta)I3 0;0 eta*I3][I3;I3]\varphi
        %   varphi = [1;0;0]+[-1 -1;1 0;0 1]*[xi1;xi2]
        %   Goal: From x compute s
        eta = (p_sample_cone(3,i)-P(3,1))/...
              (P(3,4)-P(3,1)); %z-value interpolated between extremes
        M = P*...
            [(1-eta)*eye(3) zeros(3,3);...
            zeros(3,3)  eta*eye(3)]*...
            [eye(3);eye(3)];
        %Holds x = M*varphi = Me_{1} + M*[-1 -1;1 0;0 1]*[xi1;xi2]
        x = p_sample_cone(:,i);
        x = x - M*[1;0;0];
        M = M*[-1 -1;1 0;0 1];
        M_LI = (M'*M)^(-1)*M';
        xi   = M_LI*x; clear x;
        Lag  = [xi;eta];
        V_cone(i) = nodals*...
                [(1-Lag(3))*eye(3) zeros(3,3);...
                   zeros(3,3) Lag(3)*eye(3)]*...
                [eye(3);eye(3)]*...
                [1-Lag(1)-Lag(2);Lag(1);Lag(2)];
    
        %Update the waitbar
        waitbar(i/size(V_cone,2),wbar)
    end
    close(wbar)
    
    bbest_sl = cell(1,size(p_sample_sl,2));
    %Compute values at the sliver nodes
    wbar = waitbar(0,'Progress Through Sliver Samples');
    for i=1:size(p_sample_sl,2)
        %for j=1:size(p_sample_sl{i},2)
        bbest_sl{i} = zeros(1,size(V_sl{i},2));
        %[~,~,~,pt_min,pt_max] = check_3dfindpr(p_3d,pr,tri,p_cross,...
                                           %p_sample_cone(:,1),.01);
        %Rotate sliver samples back to their [-pi/2,pi/2] position
        %but preserve the 4th row gdl index
        p_sample_sl{i} = [R zeros(3,1); zeros(1,3) 1]*p_sample_sl{i};
        
        %Shift the sliver from starting at the cone rim to starting at the
        %R_t cylinder rim
        
        E_rho = [p_sample_sl{i}(1,:)./( p_sample_sl{i}(1,:).^2 + p_sample_sl{i}(2,:).^2 ).^(1/2);...
                 p_sample_sl{i}(2,:)./( p_sample_sl{i}(1,:).^2 + p_sample_sl{i}(2,:).^2 ).^(1/2);...
                 zeros(2,size(p_sample_sl{i}(1,:),2))];
        p_sample_sl{i} = p_sample_sl{i}+... diff is (1-lambda(z))*R_t*erho
                    [R_t*(1-R_b/R_t + (R_b-R_t)/(H*R_t)*p_sample_sl{i}(3,:)).*...
                     E_rho(1,:);...
                     R_t*(1-R_b/R_t + (R_b-R_t)/(H*R_t)*p_sample_sl{i}(3,:)).*...
                     E_rho(2,:);...
                     E_rho(3:4,:)];
        
        for j=1:size(V_sl{i},2)
        %Find prism containing point
        [prfound,bbest_sl{i}(j),~,~,~]=...
            check_3dfindpr(p_3d,pr,tri,p_cross,...
                           p_sample_sl{i}(:,j),.01,pt_min,pt_max);
        %Extract nodes of that prism and value of spline at those nodes
        nodes = pr(1:6,prfound);
        nodals= Sp(nodes);
        %Evaluate the spline at the sample point
        P = p_3d(:,nodes);
        %Check that pr is formatted as expected
        test1 = unique(P(3,1:3));
        test2 = unique(P(3,4:6));
        if (size(test1,2)~=1)||(size(test2,2)~=1)
            error([num2str(prfound) '^th prism is not ordered with its base ' ... 
                   'the first 3 entries and cap face the last 3 entries '...
                   'pr list.']);
        end
        %   Have the system
        %   x = P[(1-eta)I3 0;0 eta*I3][I3;I3]\varphi
        %   s = nodals*[(1-eta)I3 0;0 eta*I3][I3;I3]\varphi
        %   varphi = [1;0;0]+[-1 -1;1 0;0 1]*[xi1;xi2]
        %   Goal: From x compute s
        eta = (p_sample_sl{i}(3,j)-P(3,1))/...
              (P(3,4)-P(3,1)); %z-value interpolated between extremes
        M = P*...
            [(1-eta)*eye(3) zeros(3,3);...
            zeros(3,3)  eta*eye(3)]*...
            [eye(3);eye(3)];
        %Holds x = M*varphi = Me_{1} + M*[-1 -1;1 0;0 1]*[xi1;xi2]
        x = p_sample_sl{i}(1:3,j);
        x = x - M*[1;0;0];
        M = M*[-1 -1;1 0;0 1];
        M_LI = (M'*M)^(-1)*M';
        xi   = M_LI*x; clear x;
        Lag  = [xi;eta];
        V_sl{i}(j) = nodals*...
                [(1-Lag(3))*eye(3) zeros(3,3);...
                   zeros(3,3) Lag(3)*eye(3)]*...
                [eye(3);eye(3)]*...
                [1-Lag(1)-Lag(2);Lag(1);Lag(2)];   
        end
                    %Update the waitbar
        waitbar(i/size(p_sample_sl,2),wbar)
    end
        close(wbar)
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%V = [V_cone V_sl];
%bbest = zeros(1,size(V,2));
%[~,~,~,pt_min,pt_max] = check_3dfindpr(p_3d,pr,tri,p_cross,...
                                           %p_sample_cone(:,1),.01);
%H = waitbar(0,'Progress Through Volume Samples');
%for i=1:size(V,2)
    %Find prism containing point
%    [prfound,bbest(i),~,~,~]=...
%        check_3dfindpr(p_3d,pr,tri,p_cross,...
                       %p_sample_cone(:,1),.01,pt_min,pt_max);
    %Extract nodes of that prism and value of spline at those points
%    nodes = pr(1:6,prfound);
%    nodals= Sp(nodes);
    %Evaluate the spline at the sample point
%    P = p_3d(:,nodes);
    %Check that pr is formatted as expected
%    test1 = unique(P(3,1:3));
%    test2 = unique(P(3,4:6));
%    if (size(test1,2)~=1)||(size(test2,2)~=1)
%        error([num2str(prfound) '^th prism is not ordered with its base ' ... 
%               'the first 3 entries and cap face the last 3 entries '...
%               'pr list.']);
%    end
    %   Have the system
    %   x = P[(1-eta)I3 0;0 eta*I3][I3;I3]\varphi
    %   s = nodals*[(1-eta)I3 0;0 eta*I3][I3;I3]\varphi
    %   varphi = [1;0;0]+[-1 -1;1 0;0 1]*[xi1;xi2]
    %   Goal: From x compute s
%    eta = (p_sample(3,i)-P(3,1))/...
%          (P(3,4)-P(3,1)); %z-value interpolated between extremes
%    M = P*...
%        [(1-eta)*eye(3) zeros(3,3);...
%         zeros(3,3)  eta*eye(3)]*...
%        [eye(3);eye(3)];
    %Holds x = M*varphi = Me_{1} + M*[-1 -1;1 0;0 1]*[xi1;xi2]
%    x = p_sample(:,i);
%    x = x - M*[1;0;0];
%    M = M*[-1 -1;1 0;0 1];
%    M_LI = (M'*M)^(-1)*M';
%    xi   = M_LI*x; clear x;
%    Lag  = [xi;eta];
%    V(i) = nodals*...
%            [(1-Lag(3))*eye(3) zeros(3,3);...
%              zeros(3,3) Lag(3)*eye(3)]*...
%            [eye(3);eye(3)]*...
%           [1-Lag(1)-Lag(2);Lag(1);Lag(2)];
    
    %Update the waitbar
%    waitbar(i/size(V,2),H)
%end
%close(H)    
%end


