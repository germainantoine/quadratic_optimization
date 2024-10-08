
T = 10;
rho=0.2;
q0  = 2;
qT= 2;


MAX=1000;

list=zeros(MAX);
list0=1:MAX;


for N=1:MAX
    dt = T/(N+1);
    Q = ((1/T)+((2*rho)/(dt^2)))*eye(N)-(rho/(dt^2))*(diag(ones(1,N-1),1)+diag(ones(1,N-1),1)');
    L1 = min(eig(Q));
    LN = max(eig(Q));
    r=(LN-L1)/(LN+L1);
    list(N)=r;
end

semilogx(list0,list);
title('Gamma en fonction de N');
exportfig(gcf,'D:\Users\Antoine\Downloads\gamma.png','Format','png','color','cmyk');