function [ Z_sd,z_scaling ] = find_actch
%Given geometry and activation disc, find bottom face of that chamber site.
%   Find the base and tops of chambers closest to the activation
%   sites given by Z_sd in data. This code taken from genmesh 
%   lines 419-438. The output is the left ends of 
%   .5I {[C I] ... [C I]} C .5I sequence.
%   Z_sd = [.25*H .33*H .64*H .78*H]; 

% Grab the data 
[ R_b,R_t,H,~,~,epsilon_0,nu,~,n_chambers,...
           ~,n_sd,Z_sd,~,...
           ~] = data;

disc = [1 0];                                                               %This section should make Giovanni's  
disc = repmat(disc,1,n_chambers-2);                                         %sd2sez code unnecessary?
disc = [0 disc 1 0];
disc = epsilon_0*disc;                  %Thickness of the C's

inter= [0 1];
inter= repmat(inter,1,n_chambers-2);
inter= [.5 inter 0 .5];
inter= nu*epsilon_0*inter;              %Thickness of the I's.

z_ticr= cumsum(inter+disc);             %Compute right ends along z-
                                        %axis of the .5I {[C I]} C
                                        %.5I sequence.  
z_ticl= z_ticr - (inter+disc);          %Compute left endpoint of 
                                        %intervals by subtracting
                                        %out their lengths.

z_ticr= z_ticr(2*cumsum(...             %I chambers are at odd 
        ones(1,n_chambers))-1);         %positions and are 
                                        %n_chambers many of them.
z_ticl= z_ticl(2*cumsum(...             %Just grabbed their left
          ones(1,n_chambers))-1);       %endpoint.
      
                                        
ram  = [zeros(1,2*n_chambers);...               %We'll join
        repmat(cumsum(ones(1,n_chambers)),1,2)];%z_ticl/r-Z_sd
found= zeros(1,n_sd);                           %in 1st row and
                                                %record in 2nd
                                                %what chamber is.
    
for i=1:n_sd                            %Find the chamber with 
    ram(1,:) = [abs(z_ticr-Z_sd(i))...  %closest face by 
                abs(z_ticl-Z_sd(i))];   %computing distance from
                                        %both faces of chamber
    ram      = sortrows(ram');          %by magnitude of distance
    ram      = ram';                    %and sorting by that value.
                                        %The other row says which
                                        %chamber have.
    
    found(i) = ram(2,1);                %Found is the index of the
                                        %chamber closest.
                                        
                                        %Reset for the next cycle.
    ram(2,:) = repmat(cumsum(ones(1,n_chambers)),1,2);
                                        
end

Z_sd = z_ticl(found);                   %By default we will place 
                                        %the isomerization on the
                                        %bottom face of the found
                                        %interdiscal chamber.

Lambda  = @(z) R_b/R_t - (R_b-R_t)/(H*R_t).*z;
z_scaling= arrayfun(Lambda,Z_sd);       %Compute cone scaling at each 
                                        %special disk height.

end

