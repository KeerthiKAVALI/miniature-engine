function mPOD = pod_complet(uFEM,x,t,k,rho,cp,u0,qL,qR,s,theta,newpod,nsnap,threshold)

% USE PRECOMPUTED REDUCED BASIS OR COMPUTE NEW ONE
if (newpod==0 && exist('POD_RB.mat','file')==2)
    load('POD_RB')
else
    % TAKE SNAPSHOTS IN TIME
    if nsnap > length(t)
        error('Number of snapshots (%d) must not be greater than total time points (%d)',nsnap,length(t))
    end
    mask = fix(linspace(1,length(t),nsnap));
    uSNP = uFEM(:,mask);
    % SINGULAR VALUE DECOMPOSITION
    [SingVec,SingVal,~] = svd(uSNP);
    SingVal = diag(SingVal);
    weight = SingVal/sum(SingVal);
    % REDUCED BASIS. ONLY EIGVECTORS FROM EIGENVALUES GREATER THAN THRESHOLD
    RBx = SingVec(:,weight>=threshold);
    save('POD_RB','RBx')
end

% RESHAPE TO ENSURE THAT X & T ARE COLUMN VECTORS
x = x(:);
t = t(:);

% SOME VARIABLES FOR CONVENIENCE
nnod = length(x);
ntime = length(t);
dt = t(2)-t(1);

% PROJECT THE FULL ORDER PROBLEM ONTO THE REDUCED BASIS
% full order matrices
M = FEM_mat_1D(x,0,0);
C = rho*cp*FEM_mat_1D(x,0,0);
K = k*FEM_mat_1D(x,1,1);
% projection
Cred = RBx'*C*RBx;
Kred = RBx'*K*RBx;
%
sred = RBx'*M*s;
qred = RBx'*[qL zeros(ntime,nnod-2) qR]';
fred = sred + qred;
%
Gred = Cred + (theta*dt*Kred);
%
fnew = fred(:,1);

% SOLVE THE REDUCED ORDER SYSTEM OF ODE
RBt = zeros(size(RBx,2),ntime);
RBt(:,1) = RBx'*u0;
for i=2:ntime
    fold = fnew;
    fnew = fred(:,i);
    f = fold + theta*(fnew-fold);
    bred = (Cred-(1-theta)*dt*Kred)*RBt(:,i-1)+ dt*f;
    Unew = Gred\bred;
    RBt(:,i) = Unew;
end

% POD MODES
mPOD{1} = RBx;
mPOD{2} = RBt';
