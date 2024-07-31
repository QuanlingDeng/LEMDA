%% Simulation evolution plots

vidfile = VideoWriter('data.mp4','MPEG-4');
open(vidfile);
fig=figure;
%set(gca,'fontsize',24)
nsp = 200;
for i = nsp:nsp:N
    
    hold on

    tem = npsol(:,:,i-1); 
    subplot(2,3,4);axis square;hold on; imagesc(tem); caxis([min(tem(:)),max(tem(:))]); colorbar; title('n0','fontsize',14)
    tem = npusol(:,:,i) - npusol(:,:,i-1); 
    subplot(2,3,5);axis square;hold on; imagesc(tem); caxis([min(tem(:)),max(tem(:))]); colorbar; title('dmvx','fontsize',14)
    tem = npvsol(:,:,i) - npvsol(:,:,i-1); 
    subplot(2,3,6);axis square;hold on; imagesc(tem); caxis([min(tem(:)),max(tem(:))]); colorbar; title('dmvy','fontsize',14)

    tem = npsol(:,:,i); 
    subplot(2,3,1);axis square;hold on; imagesc(tem); caxis([min(tem(:)),max(tem(:))]); colorbar; title(['n, T =', num2str(i*dt)],'fontsize',14)
    tem = npusol(:,:,i)./npsol(:,:,i); 
    subplot(2,3,2);axis square;hold on; imagesc(tem); caxis([min(tem(:)),max(tem(:))]); colorbar; title('mvx','fontsize',14)
    tem = npvsol(:,:,i)./npsol(:,:,i); 
    subplot(2,3,3);axis square;hold on; imagesc(tem); caxis([min(tem(:)),max(tem(:))]); colorbar; title('mvy','fontsize',14)

    hold off
    %set(gca,'fontsize',24)
    F(i) = getframe(gcf);
    writeVideo(vidfile,F(i));
    %saveas(fig,['./1022/swfloe100/' num2str(i,'%03.f') '.png'],'png');
end
close(vidfile)
toc
