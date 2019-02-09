function [ INDEX ] = indexB( N,index )
%Convert tp basis index in single or multinomial form to the other
%   N is a row vector with each entry giving the dimension of the S^r_d
%   space used in that variable. index is either a 1x1 or a 1 x size(N,2)
%   vector with natural number entries.  If it is 1x1, indexB gives the
%   multinomial index which under lexicographic order is the index^th
%   multinomial. Similarly, if it is a multinomial, indexB gives its linear
%   position with respect to lex. Key idea of script is that under lex,
%   [w < a] for 'a' given is a disjoint union(i) of words that agree with
%   'a' exactly in first (i-1) entries and wi < ai. 
%   RECALL decision trees are ordered by lex!

n = size(N,2);
if (size(index,2)~=1)&&(size(index,2)~=size(N,2))
    disp(['index must either be 1x1 or 1x' num2str(n)]);
end

switch size(index,2)
    case 1 %Convert into a multinomial
        %   Moving from left entry to right, find the maximum
        %   entry values which keep you less than or equal to
        %   linear position index.
        
        index = index - 1;      %To put start multinomial @ index = 1;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N = [N 1];              %pad for benefit of indexing below
        INDEX = zeros(1,n);
        ram = index;
        for i=1:n
            under = prod(N(1,i+1:n+1));    %Each step in i^th index
                                           %passes this many elements
                                           %under it in the tree.
 
            INDEX(i) = floor(ram/under)+1; %Find max # passes that keeps
                                           %you under. But rmbr INDEX(i)-1
                                           %is actual # of passes bc we
                                           %start at 1.
                                           
            ram = ram - (INDEX(i)-1)*under;%Redefine ram to reflect this 
                                           %last pass.
        end
        
    case n %Convert from a multinomial
        N = [N 1]; %Pad for purpose of indexing
        INDEX = 1; %So counter starts at 1                                  %%%%%%%%%%%%%%%%%%%%
        for i=1:n
            under = prod(N(1,i+1:n+1));        %The i^th counter passes     %Loop computes number of%
                                               %this many terms underneath  %multinomials strictly less% 
                                               %it.                         %than given. INDEX = 1 just adds%
                                                                            %1 so is # <=. Two cases in code%
            INDEX = (index(i)-1)*under + INDEX;%Rmbr # completed cycles is  %are inverses of each other%
                                               %index(i)-1 bc we start at
                                               %1's in each entry.
        end
        
end

end

