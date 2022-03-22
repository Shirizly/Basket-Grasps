function plotsegment(filename)

    load(filename,'BGS','BGSn','BGDepth2','BGdepthlist','BGDepth','SSDepth');
    % figure;
    hold on;
    correc = 0;
    plot(correc+BGS/n,BGDepth2,'g','lineWidth',3,'DisplayName','Local escape');
    plot(correc+BGSn/n,BGdepthlist,'+k','MarkerSize',8,'lineWidth',3,'DisplayName','Graph search');
    plot(correc+BGS/n,BGDepth,'r--','lineWidth',2,'DisplayName','Double-support');
    plot(correc+BGS/n,SSDepth,'b--','lineWidth',2,'DisplayName','Single-support');
    title('Basket Depth for Basket grasp segment - Convex polygon','interpreter','latex','fontSize',14);
    xlabel('S','interpreter','latex');
    ylabel('$\Delta U$','interpreter','latex');
    legend()
    hold off

end