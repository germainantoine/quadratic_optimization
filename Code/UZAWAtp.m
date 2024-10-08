function [X,lambda,I, RESI_lambda, RESI_X] = UZAWAtp(Q,C,B,f,lam0,TAU,ETA,EPS,IMAX)

% Methode de Uzawa.
%
% Arguments de la fonction UZAWA:
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
lambda    = lam0;
X=2*ones(length(B),1);
oldX=X;

I    = 0;

RESI_lambda = [];
RESI_X = [];


%Calcul de X0
[X,~,] = GC(Q,B-C'*lambda,oldX,ETA,IMAX);
NX=norm(X-oldX);


%Debut de la boucle
while norm(max(0,C*X-f)) >= ETA & NX >= EPS
    oldX=X;
    oldlambda=lambda;
    lambda=max(lambda+TAU*(C*X-f),0);
    [X,~,] = GC(Q,B-C'*lambda,oldX,ETA,IMAX);
    NX = norm(X-oldX);
    Nlamb=norm(lambda-oldlambda);

    
    if NX > 1E10
       warning('Explosion dans UZAWA!');
    end
    
    
    RESI_X=[RESI_X;NX];
    RESI_lambda=[RESI_lambda;Nlamb];
    I=I+1;
end
if I == IMAX
   warning('IMAX a ete atteint dans UZAWA. Augmenter IMAX!!!');
end

end