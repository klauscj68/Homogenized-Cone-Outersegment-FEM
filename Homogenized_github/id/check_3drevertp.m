function [ p ] = check_3drevertp( p )
%Map the cone nodal points back to a uniform R_t cylinder geometry
%   We assume that p and pr were generated from the diffusion codes which
%   means initially they live in the physical cone, while the samples were
%   generated in a R_t cylinder with a + epsilon_0 extension over half the
%   lip circumference. Because the code is simpler in the R_t domain, we
%   map the p back to that domain. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the data
R_b = 3;
R_t = 1;
H   = 15;
nu  = 1;
n_chambers = 100;
epsilon_0 = H/((1+nu)*(n_chambers-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_p  = size(p,2);

%Map p back into the R_t and (in sliver) R_t + eps0*e_{rho} spaces
lambda = @(z) R_b/R_t - (R_b-R_t)/(H*R_t)*z;
for i=1:n_p
    if p(1,i)^2+p(2,i)^2<=R_t^2*lambda(p(3,i))^2
        %   In the cone
        p(1:2,i) = 1/lambda(p(3,i))*p(1:2,i);
    else
        %   In the sliver
        e_rho = p(1:2,i)/sqrt(sum(p(1:2,i).^2));
        p(1:2,i) = p(1:2,i)+(1-lambda(p(3,i)))*R_t*e_rho;
    end
end


end

