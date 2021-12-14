%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MoRed - ICI/HPC Institute
% ECOLE CENTRALE DE NANTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

addpath(genpath('shared'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ENTRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same problem as Lab1_1 but this time the thermal conductivity is a
% nonlinear function of the temperature T 

% FIXED POINT TOLERANCE
fp_tol = 1e-7;

% LOAD FINITE ELEMENT MESH
load('Mesh_M2')
X = nodesAll(:,1); Y = nodesAll(:,2);
numOfNodes = size(nodesAll,1);
load('Reduced_Basis.mat')

% TIME INTERVAL, TIME INCREMENT, 
timeSpan = 1; timeInc = 0.01;
time = linspace(0,timeSpan,timeSpan/timeInc+1)'; numOfFrames = length(time);

% INITIAL CONDITION 
solution = zeros(size(reducedBasis,2),numOfFrames);

% NEWMANN BC
nodesNeu = unique(connectHole(:));
ampliNeu = sin(2*pi/timeSpan*time);
ampliNeu(ampliNeu<0) = 0;
valuesNeu = 1*ones(size(nodesNeu))*ampliNeu.';

% HEAT SOURCE
ampliSource = sin(2*pi/timeSpan*(time-timeSpan/2));
ampliSource(ampliSource<0) = 0;
valuesSource = 3*ones(numOfNodes,1)*ampliSource.';

% SOLVE: LINEAR FINITE ELEMENTS & SINGLE STEP 1st ORDER TIME INTEGRATION

% finite element operators
opMass = FEM_mat_2D(nodesAll,connectAll,0,0,0,0);
opMassNeuBc = FEM_mat_1D(nodesAll,connectHole,0,0);



% right-hand side
auxSp = allcomb(1:numOfFrames,nodesNeu);
auxNeumann = sparse(auxSp(:,2),auxSp(:,1),valuesNeu(:),numOfNodes,numOfFrames);

rhsNeumannFree = opMassNeuBc*auxNeumann;
rhsSourceFree = opMass*valuesSource;

% sum contributions and project on reduced basis
rhsTotalFree = reducedBasis'*(rhsNeumannFree + rhsSourceFree);
% Auxiliary variable used to hold Temperature fields
auxT = solution(:,1);
% TIME INTEGRATION LOOP (first order single step)
tic,
for i = 2:numOfFrames
    % 
    it = 0;
    err = 1;  
    rhs = density*specificHeat*M_r*solution(:,i-1) + timeInc*rhsTotalFree(:,i);
    while (err > fp_tol)&&(it<100)
        it = it +1;
        % evaluate the conductivity field 
        T_eim = Int_t*auxT;    % evaluate temperature at eim points
        cond_eim = condfun(T_eim); % evaluate conductivity from T_eim
        beta_cond =Int_c/cond_eim ; % evaluate beta coefficients for conductivity
        % assemble reduced LHS operator
        A = density*specificHeat*M_r; % mass matrix
        for j = 1:numel(K_r) % add other terms of affine decomposition
            A = A + timeInc* K_r{j}*beta_cond(j) ;
        end
        % solve system
        solution(:,i) = A\rhs;
        err = norm(solution(:,i)-auxT)/norm(solution(:,i));
        fprintf('Time step %d of %d , error est : %g\n',i,numOfFrames,err);
        auxT=solution(:,i);
    end
end
execTime = toc;

fprintf('Time integration loop executed in %.4f sec\n',execTime)

reducedSolution = solution;
solution = reducedBasis*solution;
save('ROM_Solution','solution','reducedSolution');

% PLOT RESULTS
figure, subplot(1,2,1), plot(time,ampliNeu,'.-')
        axis square, grid on, xlabel('time'), ylabel('neumann term amplitude'),
        set(gca,'FontSize',16)
        %
        subplot(1,2,2), plot(time,ampliSource,'.-')
        axis square, grid on, xlabel('time'), ylabel('source term amplitude'),
        set(gca,'FontSize',16)

% PLOT TIME SOLUTION
figure,
solution2plot = solution(:,1);
f1 = patch(X(connectAll.'),Y(connectAll.'),solution2plot(connectAll.'),'EdgeColor',[0.8 0.8 0.8]);
    axis equal, axis off, set(gca,'FontSize',16),
    title(['time = ' num2str(time(1),'%.2f')])
    colorbar, caxis([min(solution(:)) max(solution(:))])
for i = 1:numOfFrames
    solution2plot = solution(:,i);
    set(f1,'CData',solution2plot(connectAll.'))
    title(['time = ' num2str(time(i),'%.2f')])
    drawnow
    pause(5/numOfFrames)
end