function [ tree ] = lexB( N )
%List all multinomial indices in lex order with index N
%   
n = size(N,2);
tree = zeros(prod(N),n);

for i=1:prod(N)
    tree(i,:) = indexB(N,i);
end
end

