
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MoRed
% ECOLE CENTRALE DE NANTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ENTRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAR LENGTH AND TIME INTERVAL
L = 1; T = 40;

% MATERIAL PROPERTIES (diffusivity, mass density, specific heat)
k = 1; rho = 1; cp = 50; 

% ELEMENT DIMENSION, TIME INCREMENT
dx = 0.02; dt = 0.08; 

% TIME INTEGRATION PARAMETER FROM 0 (FULLY EXPLICIT) TO 1 (FULLY IMPLICIT)
theta = 1;

% INITIAL CONDITION 
f_u0 = @(aux) zeros(length(aux),1);
%f_u0 = @(aux) aux(:);

% DIRICHLET BC (L: left, R: right)
f_uL = @(aux) [];
f_uR = @(aux) [];

% NEWMANN BC (L: left, R: right)
f_qL = @(aux) zeros(length(aux),1);
f_qR = @(aux) ones(length(aux),1);

% HEAT SOURCE
f_s = @(auxx,auxt) zeros(length(auxx),length(auxt));
%f_s = @(auxx,auxt) ones(length(auxx),length(auxt));


%%%%%%%%%%%%%%% PREVIOUS COMPUTATIONS AND VERIFICATIONS %%%%%%%%%%%%%%%%%%%

x = linspace(0,L,L/dx+1)';
t = linspace(0,T,T/dt+1)';
u0 = f_u0(x);
uL = f_uL(t);
uR = f_uR(t);
qL = f_qL(t);
qR = f_qR(t);
s = f_s(x,t);

for i=1:length(t)
    if t(i)<=10
        qR(i)=t(i);
    elseif t(i)<=20
        qR(i)=10;
    elseif t(i)<=30
        qR(i)=30-t(i);
    else
        qR(i)=0;
    end
end

% SOME VERIFICATIONS
nnod = length(x);
ntime = length(t);
if length(u0)~=nnod
    error('Bad definition of function f_u0')
elseif size(s,1)~=nnod
    error('Bad definition of function f_s')
elseif length(uL)~=ntime && ~isempty(uL)
    error('Bad definition of function f_uL')
elseif length(uR)~=ntime && ~isempty(uR)
    error('Bad definition of function f_uR')
elseif length(qL)~=ntime
    error('Bad definition of function f_qL')
elseif length(qR)~=ntime
    error('Bad definition of function f_qR')
end


%%%%%%%%%%%%%%%%%%%%%%%%%% POD REDUCED MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET 1 TO COMPUTE NEW POD BASIS. SET 0 TO USE PREVIOUS ONE
newpod = 1;
er=[];
nsnap =[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ];
threshold = [10^-1 10^-2 10^-3 10^-4 10^-5 10^-6 10^-7 10^-8 10^-9 10^-10 10^-11 10^-12 10^-13 10^-14 10^-15 10^-16 10^-17 10^-18 10^-19 10^-20];

for n =1:20
    for m =1:20
    % POD COMPUTATION
    load('Full_Solution.mat')
    mPOD = pod_complet(uFEM,x,t,k,rho,cp,u0,qL,qR,s,theta,newpod,nsnap(n),threshold(m));
    a = norm(uFEM'-(mPOD{2}*mPOD{1}'),2);
    b = norm(uFEM',2);
    er1 = a/b ; 
    er(n,m) = er1;
    end
end
figure(1)
plot(nsnap,er)
set(gca,'Xscale','log','Yscale','log')
xlabel('log(snapshots)')
ylabel('log(error)')
title("Snapshots vs Relative error")
grid on

figure(2)
plot(threshold,er)
set(gca,'Xscale','log','Yscale','log')
xlabel('log(threshold)')
ylabel('log(error)')
title("Threshold vs Relative error")
grid on

figure(3)
[S,T] = meshgrid(nsnap,threshold);
surf(S,T,er);
xlabel('log(snapshots)')
zlabel('log(error)')
ylabel('log(threshold)')
set(gca,'Xscale','log','Zscale','log','Yscale','log')
title("Relative Error")
grid ON
