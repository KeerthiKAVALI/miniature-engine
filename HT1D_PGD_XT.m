% PGD solution of the transient heat equation with homogeneous Dirichlet boundary conditions
% and homogeneous initial condition 

function [alpha,B0s,B1s] = HT1D_PGD_XT(x,t,k,rho,cp,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,f0,f1)

% Mesh definition for each dimension
Nt = numel(t);
Nx = numel(x);
% Reshape to ensure that x & t are column vectors
t = t(:);
x = x(:);
% Mesh size
dt = (t(2)-t(1));
hx = (x(2)-x(1));

% Finite difference matrices in x
Is = rho*cp*speye(Nx);
% Finite difference second order differentiation matrix
K = k*gallery('tridiag',Nx,(1/hx)^2, -2*(1/hx)^2, (1/hx)^2);


% Finite difference matrices in time
It = k*speye(Nt);
% Time differentiation matrix
G = rho*cp*gallery('tridiag',Nt,-1/dt, 1/dt, 0);

% Containers to store the B1s and B2s functions
alpha = [];
B0s = [];
B1s = [];


% Enrichment loop
for term = 1:Max_terms
    % Initialization of the fixed point loop
    
    B0 = randn(Nx,1);
    B1 = randn(Nt,1);
    rng(1);
    
    % Satisfaction of the homogeneous Dirichlet Boundary conditions for the
    % enrichments
    B0(1,1) = 0;
    B0(end,1) = 0;
    B1(1,1) = 0;
    
    % Fixed point iterations
    for iter = 1:Max_fp_iter
        % Store the old values of B0 & B1 for later comparison
        B0_Old = B0;
        B1_Old = B1;
          
        % Solve for B0
        % Construction of the boundary value problem along x
        % LHS coefficients
        alpha_0  = trapz(t,B1.*((G./(rho*cp))*B1));
        beta_0 = trapz(t,B1.*B1);
        
        % Construction of the FD boundary value problem
        A0 = alpha_0*Is - beta_0*K;
        
        % Source term coefficient
        gamma_0 = trapz(t,f1.*B1);
        % Construction of the RHS
        b0 = f0*gamma_0;
        % In case this is not the first enrichment, previous terms are added
        % to the RHS
        if (term>1)
            % RHS coefficients
            alpha_0_i = trapz(t,bsxfun(@times,B1,(G./(rho*cp))*B1s));
            beta_0_i = trapz(t,bsxfun(@times,B1,B1s));
            b0 = b0 - (Is*(alpha.*B0s))*alpha_0_i' + (K*(alpha.*B0s))*beta_0_i';
        end
        
        % Solution with homogeneous boundary conditions
        B0(2:end-1) = A0(2:end-1,2:end-1)\b0(2:end-1);

        % Solve for B1
        % Construction of the initial value problem along t
        % LHS coefficients
        alpha_1 = trapz(x,B0.*B0);
        beta_1  = trapz(x,((K./k)*B0).*B0);
        
        % Construction of the initial value problem
        A1 = alpha_1*G - beta_1*It;
        
        % Source term coefficient
        gamma_1 = trapz(x,f0.*B0);

        % Construction of the RHS
        b1 = f1*gamma_1;
        % In case this is not the first enrichment, previous terms are added to the RHS
        if (term>1)
            % RHS coefficients
            alpha_1_i = trapz(x,bsxfun(@times,B0,B0s));
            beta_1_i = trapz(x,bsxfun(@times,B0,(K./k)*B0s));
            b1 = b1 - (G*(alpha.*B1s))*alpha_1_i' + (It*(alpha.*B1s))*beta_1_i';
        end
        
        % Solution with homogeneous initial condition.
        B1(2:end) = A1(2:end,2:end)\b1(2:end);
               
        % Norm of the difference between the 2 fixed point iterations
        S_difference = sqrt(trapz(x,B0.^2)*trapz(t,B1.*B1) + trapz(x,B0_Old.^2)*trapz(t,B1_Old.^2) - 2*trapz(x,B0.*B0_Old)*trapz(t,B1.*B1_Old));
        % Fixed point exit test
        if(S_difference < epsilon), break; end
    end
    % Normalized modes are added to the  existing ones 
    fact_0 = sqrt(trapz(x,B0.^2));
    fact_1 = sqrt(trapz(t,B1.^2));
    alpha = fact_0*fact_1;
    B0s = [B0s B0/fact_0];
    B1s = [B1s B1/fact_1];
    
    % Simplified stopping criterion
    E = alpha() / alpha_0;
    
    if(E<epsilon_tilde), break; end
    
end
end
