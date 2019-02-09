function [x,y]=incisure(bs,s) 
% disc with sliver

% legge i dati sulla geometria
[R_b,~, ~, ~, ~, ~, ...
    theta_in,theta_fin,...    
    ~, ~]=data;



% number of arcs of the circular domain contour
nbs=2;



% pass nbs
if nargin==0  
  x=nbs;   
  return 
end

% Matrix dl
dl=zeros(4,nbs);

if theta_in~=0
    disp('theta_in must be equal to 0')
end

% The arcs are parametrized by means of their initial and final angle
dl(1:2,1:nbs)=[[theta_in theta_fin]; [theta_fin 2*pi+theta_in]];


% labels assigments for zones delimited by the arcs: 1 for inner zone and 0
% for outer zone

% lati del poligono centrale = a sin zona 1, a dx zona i+1
dl(3:4,1:nbs)=[ones(1,nbs); zeros(1,nbs)];



% pass dl
if nargin==1   
  x=dl(:,bs);   
  return 
end 

% nodal positions x(s),y(s) requested in s
[m,n]=size(s);
x=zeros(m,n);
y=zeros(m,n);
[mm,nn]=size(bs); 
if mm==1 && nn==1,   
  bs=bs*ones(m,n); % expand bs 
elseif mm~=size(s,1) || nn~=size(s,2),   
  error('bs must be scalar or of same size as s'); 
end 

for ii=1:m
    for jj=1:n
        
            % archi esterni
            x(ii,jj)=R_b*cos(s(ii,jj));
            y(ii,jj)=R_b*sin(s(ii,jj));
        
    end
end

return
