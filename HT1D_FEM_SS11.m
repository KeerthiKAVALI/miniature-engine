function U = HT1D_FEM_SS11(x,t,k,rho,cp,u0,uL,uR,qL,qR,s,theta)

nnod = length(x);
ntime = length(t);
dx = x(2)-x(1);
dt = t(2)-t(1);

U = zeros(nnod,ntime);

if theta < 0.5
    dt_cr = 2/(1-2*theta)/3/k*rho*cp*dx^2;
    if  dt > dt_cr
        warning('Time increment (%f) should be less than critical value (%f)',dt,dt_cr);
    end
end

M = FEM_mat_1D(x,0,0);
C = rho*cp*FEM_mat_1D(x,0,0);
K = k*FEM_mat_1D(x,1,1);

G = C + theta*dt*K;
fnew = M*s(:,1) + [qL(1) zeros(1,nnod-2) qR(1)]';
mask_d = true(nnod,1);
if ~isempty(uL)
    mask_d(1) = false;
    U(1,:) = uL(:)';
end
if ~isempty(uR)
    mask_d(end) = false;
    U(end,:) = uR(:)';
end

U(:,1) = u0;

for i=2:ntime
    fold = fnew;
    fnew = M*s(:,i) + [qL(i) zeros(1,nnod-2) qR(i)]';
    f = fold + theta*(fnew-fold);
    b = (C-(1-theta)*dt*K)*U(:,i-1) + dt*f;
    U(mask_d,i) = G(mask_d,mask_d)\(b(mask_d) - G(mask_d,not(mask_d))*b(not(mask_d)));
end
end
