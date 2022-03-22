function fig = drawPolyGraspMoved(PG,axle,theta,d,sf1,sf2,f1pos,f2pos)
        fig = figure;
        fig = PG.drawPolygonMoved(fig,axle,theta,d,'k');
        hold on
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        N = 0.4;
        n1 = [0;N];
        if ~isempty(sf1)
        n1 = -R*PG.get('normal',sf1);
        end
        
%         t = linspace(0,2*pi);
        cent1 = f1pos-0.5*n1;
        plot(f1pos(1),f1pos(2),'.k','markersize',10)
        text(cent1(1),cent1(2),'$O_1$','fontsize',12,'interpreter','latex')
%         plot(cent1(1)+R*cos(theta),cent1(2)+R*sin(theta),'k')
        normal1 = [f1pos, f1pos + 5*n1];
%         plota(normal1,'r--');
        n2 = [0;N];
        if ~isempty(sf2)
        n2 = -R*PG.get('normal',sf2);
        end
        cent2 = f2pos-0.5*n2;
        plot(f2pos(1),f2pos(2),'.k','markersize',10)
        text(cent2(1),cent2(2),'$O_2$','fontsize',12,'interpreter','latex')
%         plot(cent2(1)+R*cos(theta),cent2(2)+R*sin(theta),'.r')
        normal2 = [f2pos, f2pos + 5*n2];
%         plota(normal2,'r--');
    end
    