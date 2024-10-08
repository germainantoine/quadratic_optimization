% Définition des données du problème
T = 1000;      % temps final
N = 1000;    % nombre de points de la grille spatiale
dt = T/(N+1); % pas de temps
rho = 2;     % paramètre matériel
q0  = 2;     % condition initiale
qT= 2;       % condition finale

% Définition de la matrice Q et du vecteur b
Q = ((1/T)+((2*rho)/(dt^2)))*eye(N)-(rho/(dt^2))*(diag(ones(1,N-1),1)+diag(ones(1,N-1),-1));
b = -(rho/(dt^2))*[q0 ; zeros(N-2,1); qT];

% Calcul de la constante c
c = (rho/(2*(dt^2)))*(qT^2+q0^2);

% Définition de la condition initiale
x0 = 2*ones(N,1);

% Paramètres pour la méthode de Gauss-Seidel
tolerance=1e-5; % tolérance de convergence
max_iterations=N*2; % nombre maximal d'itérations

% Résolution du système linéaire avec la méthode de Gauss-Seidel
tic;
[x,iter,residuals] = GCtp(Q,-b,x0,tolerance,max_iterations);  
toc;

% Affichage des résultats
disp(strcat('Convergence en :', int2str(iter),' iterations. ', 'Residu final:',num2str(residuals(end))));

% Tracé de la solution numérique obtenue
t = 0:dt:T; % vecteur temps
q = [q0;x;qT]; % solution numérique
plot(t,q);
title('N=1000  T=1000  σ=2');
exportfig(gcf,'D:\Users\Antoine\Downloads\ex1_1000_6.png','Format','png','color','cmyk'); % export de la figure en format png