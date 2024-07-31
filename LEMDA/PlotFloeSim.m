%% Simulation evolution plots
Dim_Grid = 20;
[xx,yy] = meshgrid(linspace(0,2*pi,Dim_Grid), linspace(0,2*pi,Dim_Grid));
x_vec = [reshape(xx,[],1), reshape(yy,[],1)];

vidfile = VideoWriter('F18sw.mp4','MPEG-4');
open(vidfile);
fig=figure;
set(gca,'fontsize',24)
nsp = 0.05/dt;
for i = 1:(N/nsp)
    %     clf
    plot(0,0)
    axis square
    hold on
    for l = 1:L
        th = 0:pi/50:2*pi;
        xunit = radius(l) * cos(th) + x(l,1+nsp*(i-1));
        yunit = radius(l) * sin(th) + y(l,1+nsp*(i-1));
        if(l<=Ll)
            h = plot(xunit, yunit,'-r','LineWidth',1);
        else
            h = plot(xunit, yunit,'-b','LineWidth',1);
        end
        
        %hh = plot(xunit, yunit,'color',[.2*radius(l),0.5,0.5]);
        %text(x(l,1+nsp*(i-1)),y(l,1+nsp*(i-1)),num2str(round(omega(l,1+nsp*(i-1)),1)));
        %text(x(l,1+nsp*(i-1)),y(l,1+nsp*(i-1)),num2str(round(thickness(l),1)));
        text(x(l,1+nsp*(i-1)),y(l,1+nsp*(i-1)),num2str(l));
    end
    xlim([0, 2*pi])
    ylim([0, 2*pi])
    box on
    %title(['t = ', num2str(dt*nsp*(i-1))]);
    title(['t = ', num2str(dt*nsp*( round((i-1)/4)*4 ))]);
    
    u = exp(1i * x_vec * kk) * (uhat(:,1+nsp*(i-1)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (uhat(:,1+nsp*(i-1)) .* transpose(rk(2,:)));

    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1,'color',[0.4843 0.8657 0.99882])
    %pause(0.1);
    hold off
    set(gca,'fontsize',24)
    F(i) = getframe(gcf); 
    writeVideo(vidfile,F(i));
    %saveas(fig,['./1022/swfloe100/' num2str(i,'%03.f') '.png'],'png');
end
close(vidfile)
toc
