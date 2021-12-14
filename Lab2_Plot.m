
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MoRed - ICI/HPC Institute
% ECOLE CENTRALE DE NANTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

fem = load('FEM_Solution','X','Y','connectAll','solution','time');
rom = load('ROM_Solution','solution','reducedSolution');

X = fem.X; Y = fem.Y; connectAll = fem.connectAll; time = fem.time;

variable = abs(fem.solution - rom.solution);
error = sqrt(sum(variable.^2,1));

figure,
plot(time,error,'r.-'),
    axis square, grid on, xlabel('time'), ylabel('error')
    set(gca,'FontSize',16),

figure,
variable2plot = variable(:,1);
f1 = patch(X(connectAll.'),Y(connectAll.'),variable2plot(connectAll.'),'EdgeColor',[0.8 0.8 0.8]);
    axis equal, axis off, set(gca,'FontSize',16),
    title(['time = ' num2str(time(1),'%.2f')])
    colorbar, caxis([min(variable(:)) max(variable(:))])
for i = 1:length(time)
    variable2plot = variable(:,i);
    set(f1,'CData',variable2plot(connectAll.'))
    title(['time = ' num2str(time(i),'%.2f')])
    drawnow
    pause(5/length(time))
end
