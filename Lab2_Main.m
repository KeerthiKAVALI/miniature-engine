
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
load('Mesh_M2.mat')
X = nodesAll(:,1); Y = nodesAll(:,2);
numOfNodes = size(nodesAll,1);

% TIME INTERVAL, TIME INCREMENT, 
timeSpan = 1; timeInc = 0.01;
time = linspace(0,timeSpan,timeSpan/timeInc+1)'; numOfFrames = length(time);

% MATERIAL PROPERTIES (conductivity, mass density, specific heat)
conductivity = 1; density = 1; specificHeat = 10; 
% conductivity  handle function  
condfun = @(T)(conductivity*(1+0.5*((T+1).^5-1)));
figure
plot(linspace(0,1),condfun(linspace(0,1)))
xlabel('Temperature')
ylabel('Conductivity')
% INITIAL CONDITION 
initialSolution = zeros(numOfNodes,1);

% DIRICHLET BC
nodesDir = unique([connectRight(:); connectUp(:)]);
ampliDir = zeros(numOfFrames,1);%2+sin(pi/2/timeSpan*time);%
ampliDir(ampliDir<0) = 0;
solutionDir = ones(size(nodesDir))*ampliDir.';

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

% dirichlet mask
maskDir = true(numOfNodes,1);
maskDir(nodesDir) = false;

% solution containers
solutionDir(:,1) = initialSolution(~maskDir);
solutionFree = nan(sum(maskDir),numOfFrames);
solutionFree(:,1) = initialSolution(maskDir);

% right-hand side operators (constant throughout the time integration)
opRhs = density*specificHeat*opMass;
opRhsFree = opRhs(maskDir,maskDir);
% right-hand side
auxSp = allcomb(1:numOfFrames,nodesNeu);
auxNeumann = sparse(auxSp(:,2),auxSp(:,1),valuesNeu(:),numOfNodes,numOfFrames);

rhsNeumannFree = opMassNeuBc(maskDir,:)*auxNeumann;
rhsSourceFree = opMass(maskDir,:)*valuesSource;

% sum contributions
rhsTotalFree = rhsNeumannFree + rhsSourceFree;
% Auxiliary variable used to hold Temperature fields
auxT = nan(numOfNodes,1);
% TIME INTEGRATION LOOP (first order single step)
tic,
for i = 2:numOfFrames
    % 
    err = 1;
    auxT(~maskDir)=solutionDir(:,i);
    auxT(maskDir)=solutionFree(:,i-1);
    auxRhsFree = rhsTotalFree(:,i);
    rhs = opRhsFree*solutionFree(:,i-1) + timeInc*auxRhsFree;
    while err > fp_tol
        % evaluate the conductivity field 
        condfield = condfun(auxT);
        Kx = FEM_mat_2D(nodesAll,connectAll,1,0,1,0,condfield);
        Ky = FEM_mat_2D(nodesAll,connectAll,0,1,0,1,condfield);
        opConduc = Kx+Ky;
        opLhs = density*specificHeat*opMass + timeInc*opConduc;
        opLhsFree = opLhs(maskDir,maskDir);
        solutionFree(:,i) = opLhsFree\rhs;
        err = norm(solutionFree(:,i)-auxT(maskDir))/norm(solutionFree(:,i));
        fprintf('Time step %d of %d , error est : %g\n',i,numOfFrames,err);
         auxT(maskDir)=solutionFree(:,i);
    end
end
execTime = toc;

fprintf('Time integration loop executed in %.4f sec\n',execTime)

solution = nan(numOfNodes,numOfFrames);
solution(~maskDir,:) = solutionDir;

solution(maskDir,:) = solutionFree;
save('FEM_Solution','X','Y','connectAll','time','solution','maskDir',...
                    'density','conductivity','condfun','specificHeat');

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
