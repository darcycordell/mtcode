function [id1,id0,niter] = plot_misfit_convergence(invlog)
%Function to plot the converange of an inversion for all iterations
%
% Usage: plot_misfit_convergence(invlog)
%
% Inputs: "invlog" is a structure variable containing iteration information
% from ModEM or WSINV inversions.
%
% Outputs: 
%   id1 is the iteration number of the "best trade-off" inversion (i.e. the
%   iteration which is nearest to an rms of 1 and a model norm of 0)
%
%   id0 is the iteration number of the "smallest objective function" (i.e.
%   the iteration which is nearest to an rms of 0 and a model norm of 0.
%   This *should* ALWAYS be the last iteration in theory)
%
%   niter is the total number of iterations. In general niter == id0.
%
%%
set_figure_size(1);

niter = invlog.niter;

subplot(2,4,3)
semilogy(invlog.lambda,'k.')
title('Regularization (\lambda)')

subplot(2,4,4)
semilogy(invlog.alpha,'k.')
title('Step Size (\alpha)')

subplot(2,4,7)
semilogy(invlog.m2,'k.')
title('Model Norm (m2)')
xlabel('Iteration Number')

subplot(2,4,8)
semilogy(invlog.f,'k.')
title('Objective Function (f)')
xlabel('Iteration Number')

subplot(2,4,[1 2])
plot(invlog.rms,'k.'); hold on
title('Misfit Convergence')
xlabel('Iteration number')
ylabel('RMS misfit')
text(invlog.niter*0.75,max(invlog.rms)*0.75,['Final RMS = ',num2str(invlog.rms(invlog.niter),'%5.2f')])

%log.m2 = log.lambda.*log.m2;

subplot(2,4,[5 6])
plot(invlog.rms,invlog.m2,'.k'); hold on
xlabel('RMS Misfit')
ylabel('Model Norm')
%axis equal
if ~isnan(invlog.m2)
    axis([0 max(invlog.rms)*1.5 0 max(invlog.m2)*1.5])   
end

d0 = sqrt((invlog.rms).^2+(invlog.m2).^2);
[~,id0] = min(d0);

d1 = sqrt((invlog.rms-1).^2+(invlog.m2).^2);
[~,id1] = min(d1);

%plot([0 invlog.rms(id0)],[0 invlog.m2(id0)],'--b')
%plot([1 invlog.rms(id1)],[0 invlog.m2(id1)],'--b')
%plot(invlog.rms(id1),invlog.m2(id1),'*r')
%text(max(invlog.rms)*0.75,max(invlog.m2)*0.75,['"Best" Iteration = #',num2str(id1),'. (RMS = ',num2str(invlog.rms(id1)),')']);

title('"L Curve"')

annotation('textbox', [0 0.9 1 0.08], ...
'String', ['Total Number of Iterations: ',num2str(invlog.niter)], ...
'EdgeColor', 'none', ...
'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')
