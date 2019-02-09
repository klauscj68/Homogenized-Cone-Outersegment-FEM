function [ pairs ] = check_genpts( rinc,thetainc )
%With these increments, generate the sample points in the disc.
%   Generate x-y coordinates of sample points in unit disc which are evenly
%   spaced in rinc and thetainc with given spacing.  

%Redefine so integer many points. The n is # of intervals
nr = floor(1/rinc);
rinc = 1/nr;

ntheta = floor(2*pi/thetainc);
thetainc = 2*pi/ntheta;

%List all pairs of (rad,theta) by listing in lex order.
%Key idea is that under lex the collection a <= b
%may be written as disjoint of those words agreeing
%with b in first k-1 entries and <= at the kth entry
pairs = zeros(nr*ntheta,2);
pairs(1,:) = [1 1];
k = 1;
%This loop is the successor function under lex order
%iterated on [1,1] as starting. 
while k <= nr*ntheta  %Because \theta is 
                      %2pi periodic, there are as many theta points as
                      %intervals, whereas because \r is not periodic 
                      %there are nr+1 many r points except we won't
                      %sample at r=1.  Thus again have nr many rad
                      %pts and ntheta many theta points
    %If 2nd index can be cycled without changing counter
    if pairs(k,2) < ntheta 
        pairs(k+1,:) = [pairs(k,1) pairs(k,2)+1];
    elseif pairs(k,1) < nr
        pairs(k+1,:) = [pairs(k,1)+1 1];
    end
    k = k+1;
end

%Switch from indices on increments to actual increments
%This means counter indices should start at 0 and go up
%to final value (is one increment short of r = 1
%and theta = 2*pi).
pairs = [rinc*(pairs(:,1)-1) thetainc*(pairs(:,2)-1)];

%Delete the redundant points at the origin, ie
%radius = 0 and theta ranges.
pairs(2:ntheta,:) = [];

%Now map the (r,theta) pairs into the disc.
pairs = [pairs(:,1).*cos(pairs(:,2)) ...
         pairs(:,1).*sin(pairs(:,2))];
        
end

