function [X,I,RESI] = GPtp(A,B,F,X0,ETA,TAU,IMAX,T,dt,q0,qT)

% Methode de Gradient Projete.
%
% Arguments de la fonction GP:
%   - X    => Solution approchee donne par la methode. (SORTIE)
%   - I    => Nombre d'iterations pour converger. (SORTIE)
%   - RESI => Vecteur avec l'evolution du residu. (SORTIE)
%   - A    => Matrice symetrique definie positive. (ENTREE)
%   - B    => Terme linéaire. (ENTREE)
%   - F    => Vecteur de projection  sur K. (ENTREE)
%   - X0   => Iterant initiale. (ENTREE)
%   - ETA  => Petit parametre pour le test d'arret. (ENTREE)
%   - TAU  => Pas. (ENTREE)
%   - IMAX => Nombre max d'iterations de la methode. (ENTREE)


%Initialisation
tt =[0:dt:T];
X    = X0;
I    = 0;
RESI = [];
G=A*X-B;
Y=X-TAU*G;

%Calcul de X1
oldX=X;
X=max(Y,F);
G=A*X-B;
Y=X-TAU*G;
N0 = norm(X-oldX);
RESI=[RESI;N0];
I=I+1;


%Debut de la boucle
while N0 > ETA & I < IMAX
    oldX=X;
    X=max(Y,F);
    G=A*X-B;
    Y=X-TAU*G;
    N0 = norm(X-oldX);
    N1 = norm(G);
    if N0 > 1E10
       warning('Explosion dans GP!');
    end
    
    if (mod(I,10)==0)
        N0;
        N1;
        plot(tt, [q0;X;qT], 'r');
    end
    RESI=[RESI;N0];
    I=I+1;
end
if I == IMAX
   warning('IMAX a ete atteint dans GP. Augmenter IMAX!!!');
end

end