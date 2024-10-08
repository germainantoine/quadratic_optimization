% Données d'entrée
T = 7;              % Temps final
N = 100;            % Nombre de points dans le maillage
dt = T/(N+1);       % Discrétisation en temps

% Paramètres du modèle
rho = 0.01;         % Coefficient rho
q0  = 2;            % Condition initiale q(0)
qT= 2;              % Condition finale q(T)

% Matrice de discrétisation
Q = diag(((1/T)+((2*rho)/(dt*dt)))*ones(1,N))-diag((rho/(dt*dt))*ones(1,N-1),1)-diag((rho/(dt*dt))*ones(1,N-1),1)';

% Terme de droite
b = [-(rho/(dt*dt))*q0 ; zeros(N-2,1); -(rho/(dt*dt))*qT];

% Constante c
c = (1/2)*(rho/(dt*dt))*(qT*qT+q0*q0);

% Fonction d'entrée
d0=0.75;
P=mod(dt:dt:T-dt,1);
I = find((P>=0.5) & (P<=16/24)) ;
f = -inf*ones(N,1);
f(I)=d0;

% Discrétisation de la condition de Dirichlet
Id=eye(N);
C=-Id(I,:);
fU= -d0*ones(length(I),1);

% Paramètres pour l'algorithme Uzawa
tau=L1;
eta=1.e-3;          % Résidu désiré
eps=1e-5;
imaxU=5*N;
lam0=zeros(size(fU));

% Exécution de l'algorithme Uzawa
timer();
tic;
[xU,lambdaU,nU, resU_Lambda, resU_x] = UZAWAtp(Q,C,-b,fU,lam0,tau,eta,eps,imaxU);
t=toc;

% Affichage des résultats
fprintf('temps=%5.2f (s), Iterations:%5i, ||x^k-x^{k-1}||=%10.2e\n',t,nU,resU_x(end));
disp(strcat('Convergence en :', int2str(nU),' iterations. ', 'Residu final:',num2str(resU_x(end))));
tt = [0:dt:T];

hold on
plot(tt, [q0;xU;qT], 'b', tt, [-inf;f;-inf], 'k');
title('N=100  T=7 d0=3/4 σ=0.01');
exportfig(gcf,'D:\Users\Antoine\Downloads\ex31.png','Format','png','color','cmyk');