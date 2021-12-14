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
f0=ones(size(x,1),1);
f1=ones(size(t,1),1);




[alpha,B0s,B1s] = HT1D_PGD_XT(x,t,k,rho,cp,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,f0,f1);
uPGD = (alpha.*B0s)*B1s';
[X, T]=meshgrid(x,t);
subplot(1,2,1)
surf(X,T,uPGD'), grid on, hold on
xlabel('space','fontsize',10)
ylabel('time','fontsize',10)
zlabel('temperature','fontsize',10)
legend('PGD Solution','Location', 'northoutside', 'fontsize',10)

% Exact Solution
N=200;
uEXA = ExactSolution(x,t,N);
[X, T]=meshgrid(x,t);
subplot(1,2,2)
surf(X,T,uEXA'), grid on, hold on
xlabel('space','fontsize',10)
ylabel('time','fontsize',10)
zlabel('temperature','fontsize',10)
legend('Exact Solution','Location', 'northoutside', 'fontsize',10)

% Relative Error
RE=norm(uPGD'-uEXA',2)/norm(uEXA',2);
