function [X,I,RESI]=GCtp(A,B,X0,ETA,IMAX)
%
% Methode de Gradient Conjugue.
%
% Arguments de la fonction GC:
%   - U    => Solution approchee donne par la methode. (SORTIE)
%   - I    => Nombre d'iterations pour converger. (SORTIE)
%   - RESI => Vecteur avec l'evolution du residu. (SORTIE)
%   - A    => Matrice symetrique definie positive. (ENTREE)
%   - B    => Second membre pour la resolution du systeme. (ENTREE)
%   - X0   => Iterant initiale. (ENTREE)
%   - ETA  => Petit parametre pour le test d'arret. (ENTREE)
%   - IMAX => Nombre max d'iterations de la methode. (ENTREE)
%
X    = X0;
G    = A*X-B;
D    = -G;
N0   = G'*G;
I    = 0;
RESI = [];
while N0 > ETA^2 & I < IMAX
   Z    = A*D;
   RHO  = -dot(G,D)/dot(Z,D);
   X    = X + RHO*D;
   G    = G + RHO*Z;
   N1   = dot(G,G);
   BETA = N1/N0;
   N0   = N1;
   D    = -G + BETA*D;
   I    = I+1; 
   RESI = [RESI,sqrt(N1)];
   if sqrt(N1) > 1E10
       warning('Explosion dans GC!');
   end
end
if I == IMAX
   warning('IMAX a ete atteint dans GC. Augmenter IMAX!!!');
end

end