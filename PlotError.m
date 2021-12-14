close all
clear
clc




L=1;
T=0.1;
k=1;
rho=1;
cp=1;

dx=0.02;
dt=0.002;
epsilon = 1e-10;
epsilon_tilde = 1e-5;
Max_terms=50;
Max_fp_iter=120;



x=linspace(0,L,L/dx + 1)';
t=linspace(0,T,T/dt + 1)';
N=50;
f0=ones(size(x,1),1);
f1=ones(size(t,1),1);
[alpha,B0s,B1s] = HT1D_PGD_XT(x,t,k,rho,cp,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,f0,f1);

er=[];
m=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
for d=1:20
    
    uPGD = (alpha.*B0s(d))*B1s(d)';
    uEXA = ExactSolution(x,t,N);
    RE=norm(uPGD'- uEXA',2)/norm(uEXA',2);
    er=[er,RE];
    
end

plot(m,er)
semilogy(er)
xlabel('Enrichment modes')
ylabel('Relative error')
title("Enrichment modes vs Relative error")
grid on

