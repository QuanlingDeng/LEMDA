refig = findobj(gca, 'Type', 'line');
fgx = get(refig, 'Xdata');
fgy = get(refig, 'Ydata');
fgxx = fgx{1,1};
fgyy = [fgy{1,1}; fgy{2,1}; fgy{3,1}; fgy{4,1}];

%%
figure;
subplot(4,4,16)
hold on
yyaxis left
plot(fgxx, fgyy(4,:), '-', 'linewidth',2)
plot(fgxx, fgyy(3,:), '--', 'linewidth',2)
ylabel('RMSE')

yyaxis right
plot(fgxx, fgyy(2,:), '-', 'linewidth',2)
plot(fgxx, fgyy(1,:), '--', 'linewidth',2)
ylabel('Pattern correlation')

title('Comparison of the skill scores','fontsize',18)
set(gca,'fontsize',24); set(gca,'linewidth',2)
%legend('vx','vy')
box on
xlabel('t')


%%
refig = findobj(gca, 'Type', 'line');
fgx = get(refig, 'Xdata');
fgy = get(refig, 'Ydata');
fgxx = fgx;
%fgyy = fgy;
fgyy = [fgyy; fgy];

%%
figure;
subplot(1,4,4)
hold on
yyaxis left
plot(fgxx, fgyy(1,:), '-*', 'linewidth',2)
ylabel('RMSE')

yyaxis right
plot(fgxx, fgyy(2,:), '-*', 'linewidth',2)
ylabel('Pattern correlation')

title('Comparison of the skill scores','fontsize',18)
set(gca,'fontsize',24); set(gca,'linewidth',2)
%legend('vx','vy')
box on
xlabel('t')
