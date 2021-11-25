close all
clear

load('plotPOD')

% AVAILABLE VARIABLES 'x','t','uFEM','mPOD'

% RECONSTRUCT POD SOLUTION
uPOD = mPOD{1}*mPOD{2}.';

% PLOT NON REDUCED MODEL RESULT
sample = fix(linspace(1,numel(t),10));
sample(1) = 1; sample(end) = numel(t);

for i=1:numel(sample)
    subplot(1,2,1),plot(x,uFEM(:,sample(i)),'-')
    grid on, hold on
    ylim(1.05*[min(uFEM(:)) max(uFEM(:))])
    xlabel('length (x)','fontsize',18)
    ylabel('temperature','fontsize',18)
    title(['Full order t = ' num2str(t(sample(i)),'%.2f')],'fontsize',18)
    set(gca,'fontsize',18)

    % COMPARE FULL ORDER AND POD SOLUTIONS
    subplot(1,2,2),plot(x,uPOD(:,sample(i)),'-')
    grid on, hold on
    ylim(1.05*[min(uPOD(:)) max(uPOD(:))])
    xlabel('length (x)','fontsize',18)
    ylabel('temperature','fontsize',18)
    title(['POD t = ' num2str(t(sample(i)),'%.2f')],'fontsize',18)
    set(gca,'fontsize',18)
    pause(2)
end

% PLOT SPACE-TIME SOLUTION
figure
[X,T] = meshgrid(x,t);
surf(X,T,mPOD{2}*mPOD{1}'), shading interp
set(gca,'fontsize',18)
xlabel('x','fontsize',18)
ylabel('t','fontsize',18)

% PLOT SPACE AND TIME POD MODES
figure
plot(x,mPOD{1},'-')
grid on
xlabel('length (x)','fontsize',18)
title('normalized space modes','fontsize',18)
set(gca,'fontsize',18)
figure
plot(t,bsxfun(@rdivide,mPOD{2},sqrt(diag(mPOD{2}'*mPOD{2}))'),'-')
grid on
xlabel('time (t)','fontsize',18)
title('normalized time modes','fontsize',18)
set(gca,'fontsize',18)
