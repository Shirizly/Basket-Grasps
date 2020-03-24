function fig = drawPolyGraspMoved(PG,axle,theta,d,sf1,sf2,f1pos,f2pos)
        fig = figure;
        fig = PG.drawPolygonMoved(fig,axle,theta,d,'k');
        hold on
        R = 0.4;
        n1 = R*PG.get('normal',sf1).';
        
%         t = linspace(0,2*pi);
        cent1 = f1pos-n1;
        plot(f1pos(1),f1pos(2),'.k','markersize',10)
        text(cent1(1)-0.5*R,cent1(2),'$O_1$','fontsize',12,'interpreter','latex')
        % plot(cent1(1)+R*cos(t),cent1(2)+R*sin(t),'k')
        n2 = [0;R];
        if ~isempty(sf2)
        n2 = R*PG.get('normal',sf2).';
        end
        cent2 = f2pos-n2;
        plot(f2pos(1),f2pos(2),'.k','markersize',10)
        text(cent2(1)-0.5*R,cent2(2),'$O_2$','fontsize',12,'interpreter','latex')
        % plot(cent2(1)+R*cos(t),cent2(2)+R*sin(t),'k')
    end
    