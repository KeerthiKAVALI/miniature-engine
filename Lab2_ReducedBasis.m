
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MoRed - ICI/HPC Institute
% ECOLE CENTRALE DE NANTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

addpath(genpath('shared'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ENTRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD FULL ORDER SOLUTION
load('FEM_Solution')

% SUBSPACE TRUNCATION THRESHOLD
threshold = 1e-4;

% EVALUATE NON-LINEAR TERM (CONDUCTIVITY)
condfield = condfun(solution);
solution = solution(maskDir,:);

% SINGULAR VALUE DECOMPOSITION
[leftSingularVectors,singularValues,~] = svd(solution,'econ');
maskTruncation = false(size(leftSingularVectors,2),1);

% SUBSPACE TRUNCATION
singularValues = diag(singularValues);
weights = 1 - cumsum(singularValues)/sum(singularValues);
maskTruncation(weights>=threshold) = true;

% REDUCED BASIS. ONLY EIGVECTORS FROM EIGENVALUES GREATER THAN THRESHOLD
reducedBasis = zeros(length(maskDir),sum(maskTruncation));
reducedBasis(maskDir,:) = leftSingularVectors(:,maskTruncation);

%% REDUCED BASIS FOR NONLINEAR TERM
% same as above but this time it is applied to condfield
% SINGULAR VALUE DECOMPOSITION
[leftSingularVectors,singularValues,~] = svd(condfield,'econ');
maskTruncation = false(size(leftSingularVectors,2),1);

% SUBSPACE TRUNCATION
singularValues =diag(singularValues); 
weights_c =1 - cumsum(singularValues)/sum(singularValues);
maskTruncation(weights_c>=threshold) = true;

% REDUCED BASIS. ONLY EIGVECTORS FROM EIGENVALUES GREATER THAN THRESHOLD
reducedBasis_c = leftSingularVectors(:,maskTruncation);

% Assemble reduced Mass Matrix 
% Hi-Fi Operator
M = FEM_mat_2D([X,Y],connectAll,0,0,0,0);
% Galerkin Projection of the mass matrix onto the reduced basis
M_r = reducedBasis_c'*M*reducedBasis_c;

% Assemble reduced Diffusion Operator Matrix (affine decomposition)

K_r = cell(size(reducedBasis_c,2),1);
for i = 1:size(reducedBasis_c,2)    
    Kx = FEM_mat_2D([X,Y],connectAll,1,0,1,0,reducedBasis_c(:,i));
    Ky = FEM_mat_2D([X,Y],connectAll,0,1,0,1,reducedBasis_c(:,i));
    opConduc = Kx+Ky;
    K_r{i} = reducedBasis'*opConduc*reducedBasis ;%conductivity
end


% EIM POINTS

intPoints = eim(reducedBasis_c);

% Interpolation Matrices 

Int_t = reducedBasis(intPoints,:); % temperature (from alphas to point evaluation of u)
Int_c = reducedBasis_c(intPoints,:); % conductivity (from betas to point evaluation of kappa)

save('Reduced_Basis','reducedBasis','reducedBasis_c','M_r','K_r','Int_t',...
                     'Int_c','density','conductivity','specificHeat','condfun');

% PLOT RESULTS
figure, semilogy(1:length(weights),weights,'.-')
        axis square, grid on, xlabel('subspace dimension'), ylabel('weights'),
        set(gca,'FontSize',16)
        hold on
        semilogy(1:length(weights(weights>0)),threshold*ones(length(weights(weights>0)),1),'r--')

% PLOT BASIS FUNCTIONS
figure,
function2plot = reducedBasis(:,1);
f1 = patch(X(connectAll.'),Y(connectAll.'),function2plot(connectAll.'),'EdgeColor',[0.8 0.8 0.8]);
    axis equal, axis off, set(gca,'FontSize',16),
    title('basis function #1')
    colorbar,
for i = 1:size(reducedBasis,2)
    function2plot = reducedBasis(:,i);
    set(f1,'CData',function2plot(connectAll.'))
    title(['basis function #' num2str(i)])
    drawnow
    pause(5/size(reducedBasis,2))
    %pause(1)
end


% PLOT RESULTS
figure, semilogy(1:length(weights_c),weights_c,'.-')
        axis square, grid on, xlabel('subspace dimension'), ylabel('weights (conductivity field)'),
        set(gca,'FontSize',16)
        hold on
        semilogy(1:length(weights_c(weights_c>0)),threshold*ones(length(weights_c(weights_c>0)),1),'r--')

% PLOT BASIS FUNCTIONS
figure,
function2plot = reducedBasis_c(:,1);
f1 = patch(X(connectAll.'),Y(connectAll.'),function2plot(connectAll.'),'EdgeColor',[0.8 0.8 0.8]);
    axis equal, axis off, set(gca,'FontSize',16),
    title('basis function #1')
    colorbar,
for i = 1:size(reducedBasis_c,2)
    function2plot = reducedBasis_c(:,i);
    set(f1,'CData',function2plot(connectAll.'))
    title(['basis function #' num2str(i)])
    hold on
    plot(X(intPoints(i)),Y(intPoints(i)),'or','MarkerSize',10);    
    drawnow
    pause
    %pause(1)
end
