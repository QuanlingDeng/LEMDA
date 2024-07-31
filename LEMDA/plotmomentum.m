%%
mv = zeros(2,N);

for j=1:N
    mv(1,j) = sum(sum(npusol(:,:,j)));
    mv(2,j) = sum(sum(npvsol(:,:,j)));
end

figure
hold on
plot(dt:dt:N*dt, mv(1,:), '-b', 'linewidth',2)
plot(dt:dt:N*dt, mv(2,:), '-g', 'linewidth',2)
plot(dt:dt:N*dt, sqrt(0.5*mv(2,:).^2 + 0.5*mv(1,:).^2), '--r', 'linewidth',2)
title('Total momentum','fontsize',14)
set(gca,'fontsize',16)
legend('mvx','mvy', 'total |mv|')
box on
xlabel('t')
