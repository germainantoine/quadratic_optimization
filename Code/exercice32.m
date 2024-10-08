% Données d'entrée

temps = 7;
N = 10;
pas = temps/(N+1);

rho = 0.01;
q0  = 2;
qT = 2;

% Matrice Q et vecteur b
Q = diag(((1/temps)+((2*rho)/(pas*pas)))*ones(1,N))-diag((rho/(pas*pas))*ones(1,N-1),1)-diag((rho/(pas*pas))*ones(1,N-1),1)';
b = [-(rho/(pas*pas))*q0 ; zeros(N-2,1); -(rho/(pas*pas))*qT];

% Vecteur f
f = -inf*ones(N,1);
I = find((mod(pas:pas:temps-pas,1)>=0.5) & (mod(pas:pas:temps-pas,1)<=16/24));
f(I) = 0.75;

% Contraintes
Id = eye(N);
C1 = -Id(I,:);
fU1 = -0.75*ones(length(I),1);

C0 = ones(N);
C2 = pas*triu(C0)';
fU2 = 5*ones(N,1);

C = [C1 ; C2];
fU = [fU1 ; fU2];

% Paramètres de l'algorithme Uzawa
tau = 1;
eta = 1.e-3; % Résidu désiré
eps = 1e-5;
imaxU = 5*N;
lam0 = zeros(size(fU));

% Résolution du problème d'optimisation
tic;
[x, lambda, nb_iterations, residu_lambda, residu_x] = UZAWAtp(Q, C, -b, fU, lam0, tau, eta, eps, imaxU);
temps_exec = toc;

% Affichage des résultats
fprintf('Temps d''exécution = %5.2f (s), Nombre d''itérations = %5i, ||x^k-x^{k-1}|| = %10.2e\n',temps_exec, nb_iterations, residu_x(end));
disp(strcat('Convergence en :', int2str(nb_iterations), ' iterations. ', 'Residu final:', num2str(residu_x(end))));
tt = [0:pas:temps];
plot(tt, [q0;x;qT], 'b', tt, [-inf;f;-inf], 'k');
title('N = 10, T = 7, d0 = 3/4, sigma = 0.01');
exportfig(gcf, 'D:\Users\Antoine\Downloads\ex32_10.png', 'Format', 'png', 'color', 'cmyk');