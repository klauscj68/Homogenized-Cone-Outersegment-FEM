% determina gli elementi tri(i) della mesh di triangoli definita da p,t cui p_find(:,i) appartiene
% e campiona in p_find(:,i) le corrispondenti funzioni di forma F(:,i)
% questo codice � vettorializzato sui triangoli
% nota: se i punti da piazzare sono molto numerosi, conviene usare
% find_tri_large_N, che � vettorializzato sui punti
function [tri, F, pt_min, pt_max]=find_tri(p, t, p_find, R, pt_min, pt_max)

%%%%%%%%%%%%% Script appears to work when t and p have
%%%%%%%%%%%%% have additional rows after resp the first
%%%%%%%%%%%%% 3 and 2.



% le matrici pt_min e pt_max determinano i rettangoli minimi con
% lati paralleli agli assi cartesiani, contenenti ciascun triangolo
% pt_min(:,j) � il vertice in basso a sinistra del rettangolo minimo relativo al triangolo j
% pt_max(:,j) � il vertice in alto a destra del rettangolo minimo relativo al triangolo j

% il calcolo di queste matrici � alquanto dispendioso
% per� dipende dalla mesh e non da p_find
% quindi pu� essere eseguito una volta soltanto, passandole tra i dati le
% volte successive

% p_find � una matrice con 2 colonne; p_find(:,i) contiene le coordinate
% del punto da processare

% metodo:
% ogni triangolo � pensato contenuto in un rettangolo minimo con
% lati paralleli agli assi cartesiani
% i triangoli "papabili" sono quelli il cui rettangolo minimo contiene il punto
% p_find(:,i)

% tolleranza utilizzata
tol=1000*eps;

% tolleranza aggiuntiva sul bounding box, dovuta al fatto che l'unione dei
% triangoli � strettamente contenuta nel cerchio
tol_pap=tol+max([R-max(p(1,:)), R+min(p(1,:)), R-max(p(2,:)), R+min(p(2,:))]);

% determina se deve calcolarsi pt_min e pt_max
if nargin==4,                                                               %Script generates pt_min/max if not already given as args

    % calcola pt_min e pt_max

    % le coordinate dei vertici di ciascun triangolo sono:
    % primi vertici: p(:,t(1,:))
    % secondi vertici: p(:,t(2,:))
    % terzi vertici: p(:,t(3,:));

    % determina i rettangoli minimi contenenti ciascun triangolo
    % minimo di ciascuna coordinata
    pt_min=zeros(2,size(t,2));
    pt_min(1,:)=min([p(1,t(1,:));p(1,t(2,:));p(1,t(3,:))],[],1);            %Finds the min x-val of each tri's verts
    pt_min(2,:)=min([p(2,t(1,:));p(2,t(2,:));p(2,t(3,:))],[],1);            %Finds the min y-val of each tri's verts
    % massimo di ciascuna coordinata
    pt_max=zeros(2,size(t,2));
    pt_max(1,:)=max([p(1,t(1,:));p(1,t(2,:));p(1,t(3,:))],[],1);            %Finds the max x-val of each tri's verts
    pt_max(2,:)=max([p(2,t(1,:));p(2,t(2,:));p(2,t(3,:))],[],1);            %Finds the max y-val of each tri's verts
    
end




% inizializza gli array di output
n_p_find=size(p_find,2);                                                    %Number of points to find tri's for
tri=zeros(1,n_p_find);                                                      %List to store those triangles
F=zeros(3,n_p_find);

% ciclo su tutti i punti da piazzare sulla mesh
% nota: se i punti da piazzare sono molto numerosi, conviene usare
% find_tri_large_N


for i=1:n_p_find
    
    
    
    % individua i triangoli papabili

%     % perch� il punto di misura � nel rettangolo minimo che li contiene
%     pap=find( (p_find(1,i)>pt_min(1,:)-tol_pap)&(p_find(1,i)<pt_max(1,:)+tol_pap) & ...
%               (p_find(2,i)>pt_min(2,:)-tol_pap)&(p_find(2,i)<pt_max(2,:)+tol_pap) );

    % metodo short circuit
    pap=find( (p_find(1,i)>pt_min(1,:)-tol_pap));                           %Locate tri's (by cols) whose min x-vert val is beneath p-x to find
    pap=pap( (p_find(1,i)<pt_max(1,pap)+tol_pap));                          %Within above find tri's w/max x-vert val above p-x to find
    pap=pap( (p_find(2,i)>pt_min(2,pap)-tol_pap));                          %Within above now do same with p-y to find
    pap=pap( (p_find(2,i)<pt_max(2,pap)+tol_pap) );                         %Now have tri's whose gen'd rect's contain p up to tolerance.
    
    

    % individua i triangoli che contengono il punto cercato
    % nel caso in cui ci fosse pi� di un triangolo
    % (a causa della tolleranza imposta), prende il primo trovato
    % il triangolo � orientato, quindi se (i,j,k) � una permutazione positiva,
    % (pj-pi)\cross(pk-pi))\dot verosre k � positivo
    % quindi p appartiene al triangolo se 
    % (pj-pi)\cross(p-pi)\dot versore k � positivo
    % cio�
    % \cross(p-pi)\dot [versore k \cross(pj-pi) ] � positivo
    % nota: [versore k \cross(pj-pi) ] � l'hodge di (pj-pi)
    % per le tre permutazioni positive
    % 1,2,3
    % 2,3,1
    % 3,1,2
    % primo vertice dei papabili                                            %Extract vertices of tri's whose rect's contain p
    X1=p(1,t(1,pap));
    Y1=p(2,t(1,pap));
    % secondo vertice dei papabili
    X2=p(1,t(2,pap));
    Y2=p(2,t(2,pap));
    % terzo vertice dei papabili
    X3=p(1,t(3,pap));
    Y3=p(2,t(3,pap));
    % punto cercato
    Xp=p_find(1,i);                                                         %Coordinates of point to find
    Yp=p_find(2,i);
    % ricerca
    % passa dalla numerazione locale a pap (triangoli papabili) a quella
    % globale dei triangoli                                                 %Logical indexing is used to find tri's whose vertices
                                                                            %satisfy certain constraints stip'd by point 2 be found.
    tri_temp=pap(((Y2-Y1).*(Xp-X1)-(X2-X1).*(Yp-Y1)<=tol) & ((Y3-Y2).*(Xp-X2)-(X3-X2).*(Yp-Y2)<=tol) & ((Y1-Y3).*(Xp-X3)-(X1-X3).*(Yp-Y3)<=tol) );

    
    % tratta il caso cattivo in cui tri � vuoto, 
    % a causa del fatto che l'unione dei triangoli � strettamente contenuta nel
    % cerchio
    if isempty(tri_temp)                                                    %If no triangle satisfied constraints, find one by least squares
        % prende il papabile per cui � minima la somma dei quadrati delle   %distance.
        % distanze dai tre vertici
        dist_sq=(Xp-X1).^2+(Yp-Y1).^2+(Xp-X2).^2+(Yp-Y2).^2+(Xp-X3).^2+(Yp-Y3).^2;
        [dummy,tri_temp]=min(dist_sq);
        % passa dalla numerazione locale a pap (triangoli papabili) a quella
        % globale dei triangoli
        tri_temp=pap(tri_temp);
    end

    % nel caso in cui ci fosse pi� di un triangolo cui appartiene il punto in esame, 
    % (a causa della tolleranza imposta), prende il primo trovato
    tri(i)=tri_temp(1);
    
   
    

    % indici dei nodi del triangolo trovato
    n=t(1:3,tri(i));                                                        %My edit was to restrict to rows 1:3 bc of
                                                                            %how my scripts work by index.
    % campiona le funzioni di forma nel punto p_find(:,i)                   %Samples the shape functions at the point p_find
    F(:,i)=tri_F(p(1,n),p(2,n),p_find(1,i),p_find(2,i));


% % % % %     % metodo di ricerca su tutti i triangoli (non utilizza i rettangoli contenenti)
% % % % %     % non funziona se il punto � nel cerchio ma fuori dall'unione del triangoli
% % % % %     % individua i triangoli che contengono il punto cercato, utilizzando una tolleranza tol
% % % % %     % coordinate dei vertici di ciascun triangolo
% % % % %     % primo vertice dei papabili
% % % % %     X1=p(1,t(1,:));
% % % % %     Y1=p(2,t(1,:));
% % % % %     % secondo vertice dei papabili
% % % % %     X2=p(1,t(2,:));
% % % % %     Y2=p(2,t(2,:));
% % % % %     % terzo vertice dei papabili
% % % % %     X3=p(1,t(3,:));
% % % % %     Y3=p(2,t(3,:));
% % % % %     Xp=p_find(1,i);
% % % % %     Yp=p_find(2,i);
% % % % %     tri_temp=find(((Y2-Y1).*(Xp-X1)-(X2-X1).*(Yp-Y1)<=tol) & ((Y3-Y2).*(Xp-X2)-(X3-X2).*(Yp-Y2)<=tol) & ((Y1-Y3).*(Xp-X3)-(X1-X3).*(Yp-Y3)<=tol) );
% % % % %     
% % % % %     % nel caso in cui ci fosse pi� di un triangolo cui appartiene il punto in esame, 
% % % % %     % (a causa della tolleranza imposta), prende il primo trovato
% % % % %     tri(i)=tri_temp(1);
% % % % % 
% % % % % 
% % % % % 
% % % % %     % metodo di ricerca basato su tri2grid
% % % % %     % non funziona se il punto � nel cerchio ma fuori dall'unione del triangoli
% % % % %     u=zeros(size(p,2),1); % valore fittizio
% % % % %     [uxy,tri(i),a2,a3]=tri2grid(p,t,u,p_find(1,i),p_find(2,i));
% % % % %     % 
% % % % %     % indici dei nodi del triangolo trovato
% % % % %     ind=t(:,tri(i));
% % % % %     % 
% % % % %     % campiona le funzioni di forma nel punto (Xp,Yp)
% % % % %     % questo � il carico corrispondente all'integrale della delta di Dirac in (X,Y)
% % % % %     % primo metodo: chiama tri_F
% % % % %     % coordinate nodali
% % % % %     X=p(1,ind);
% % % % %     Y=p(2,ind);
% % % % % 
% % % % %     F(:,i)=tri_F(X,Y,p_find(1,i),p_find(2,i));
% % % % %     
% % % % % %     % secondo metodo: usa a2,a3 fornite da tri2grid
% % % % % %     F(:,i)=[1-a2-a3; a2; a3];
    
end

return
