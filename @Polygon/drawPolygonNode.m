function drawPolygonNode(PG,theta,sf1,sf2,f1pos,f2pos)
        hold on
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        N = 0.6;
        n1 = [0;-N];
        if ~isempty(sf1)
            n1 = -R*PG.get('normal',sf1)*N;
        end
        if n1(1)<0
            cent1 = f1pos+n1;
        else
            cent1 = f1pos+n1;
        end
        plot(f1pos(1),f1pos(2),'.k','markersize',35)
        text(cent1(1),cent1(2),'$O_1$','fontsize',24,'interpreter','latex')
        n2 = [0;-N];
        if ~isempty(sf2)
            n2 = -R*PG.get('normal',sf2)*N;
        end
        if n2(1)<0
            cent2 = f2pos+3*n2;
        else
            cent2 = f2pos+2*n2;
        end
        plot(f2pos(1),f2pos(2),'.k','markersize',35)
        text(cent2(1),cent2(2),'$O_2$','fontsize',24,'interpreter','latex')

    end