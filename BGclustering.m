function [CC, Find] = BGclustering(M,Nb,Nbrecalc)

% Entrées :     M  =        Matrice binaire d'entrée
%               Nb =        Nombre de valeur propres calculées (500 - 10min, 1000
%                           - 30min)
%               Nb recalc = Nb de nodes recalculés de manière exacte


n = length(M);
CC = zeros(n,1);
M = Binarisation(M);

opts.tol = 1e-3;
% Calcul des valeurs et vecteurs propres
[V, D]=eigs(M,Nb,'la', opts);


V = double(V);
V = V.*V;
D = D*D*D;
% A = (D<0);
% D(A) = 0;
K = sum(V*D,2);
Deg = sum(M);

% Calcul des CC individuels
for j = 1:n
    d = Deg(j);
    if d>=2
        CC(j)=K(j)/(d*d-d);
    end
end

% Rectification des CC trouvés >1 ou <0
S = CC>1;
CC(S) = 1;
S = CC<0;
CC(S) = 0;

% Calcul exact pour les CC les plus divergents

if Nbrecalc > 0
    
[v, I]= sort(abs(V(:,1)));
i = 1;

% On cherche dans les petites valeurs dans les coefficients du 
% premier vecteur propre (mais qui sont pas trop petites)

while (v(i+1) <= v(i)*1e6) && (i < n-1) 
    i = i+1;
end

if i == n-1
    fprintf('Gap dans le vecteur propre non détecté\n')
else
    Find = I(i+1:i+1+Nbrecalc);
    for j = 1:length(Find)
        progressbar(j,length(Find));
        u = Find(j);
        CC(u) = CC_Bin_Noeud(M,u);
    end
end

end
