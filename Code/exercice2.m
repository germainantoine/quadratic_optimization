% Initialisation des données
T = 7; % Temps final
N = 100; % Nombre de pas de temps
dt = T/(N+1); % Pas de temps

rho = 0.01; % Coefficient rho
q0  = 2; % Valeur initiale de q
qT = 2; % Valeur finale de q

% Calcul de la matrice Q, du vecteur b et du scalaire c
Q = diag(((1/T)+((2*rho)/(dt*dt)))*ones(1,N))-diag((rho/(dt*dt))*ones(1,N-1),1)-diag((rho/(dt*dt))*ones(1,N-1),1)';
b = [-(rho/(dt*dt))*q0 ; zeros(N-2,1); -(rho/(dt*dt))*qT];
c = (1/2)*(rho/(dt*dt))*(qT*qT+q0*q0);

% Initialisation du vecteur x0, de la valeur seuil d0 et des paramètres de convergence
x0 = 2*ones(N,1); % Vecteur initial x0
eta = 1e-5; % Tolérance pour la convergence
imax = 1000; % Nombre maximum d'itérations

% Définition de l'ensemble K 
d0 = 0.75; % Valeur seuil d0
P = mod(dt:dt:T-dt,1); % Vecteur de temps modifié
I = find((P>=0.5) & (P<=16/24)); % Indices pour lesquels P est compris entre 0.5 et 16/24
f = -inf*ones(N,1); % Initialisation du vecteur f
f(I) = d0; % Définition de f pour les indices I

% Définition des paramètres pour la méthode GP
L1 = min(eig(Q)); % Plus petite valeur propre de Q
LN = max(eig(Q)); % Plus grande valeur propre de Q
tau = 2/(L1+LN); % Pas pour la méthode GP
r = (LN-L1)/(LN+L1); % Coefficient r
eta2 = (1-r)/r*eta; % Nouvelle tolérance pour la convergence
imaxGP = N*100; % Nombre maximum d'itérations pour la méthode GP

% Calcul de xGP avec la méthode GP
timer();
tic;
[xGP,nGP,resGP] = GPtp(Q,-b,f,x0,eta2,tau,imaxGP,T,dt,q0,qT); % Appel de la fonction GP
t = toc; % Temps d'exécution
fprintf('temps=%5.2f (s), Iterations:%5i, ||x^k-x^{k-1}||=%10.2e\n',t,nGP,resGP(end)); % Affichage des résultats

disp(strcat('Convergence en :', int2str(nGP),' iterations. ', 'Residu final:',num2str(resGP(end)))); % Affichage des résultats

tt = [0:dt:T]; % Vecteur de temps
plot(tt, [q0;xGP;qT], 'r', tt, [-inf;f;-inf], 'k');
title('N=100  T=7 d0=3/4 σ=0.01');
exportfig(gcf,'D:\Users\Antoine\Downloads\ex2.png','Format','png','color','cmyk');