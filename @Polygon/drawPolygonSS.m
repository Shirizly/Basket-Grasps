function fig = drawPolygonSS(PG)
        fig = figure;
        hold on
        % edges
        dx = mean(abs(diff(PG.vertex(1,:))));
        dy = mean(abs(diff(PG.vertex(2,:))));
        plot(PG.vertex(1,:),PG.vertex(2,:),'k')
        plot(PG.vertex(1,[1,end]),PG.vertex(2,[1,end]),'k')
%         for k=1:PG.nv-1
%             text((PG.vertex(1,k)+PG.vertex(1,k+1))/2+dx/6-1*dx*((PG.normal(1,k)>0)),(PG.vertex(2,k)+PG.vertex(2,k+1))/2-dy/10*PG.normal(2,k),['$e_{' num2str(k) '}$'],'interpreter','latex','fontsize',12)
%         end
%         text((PG.vertex(1,end)+PG.vertex(1,1))/2-dx/10*PG.normal(1,end),(PG.vertex(2,end)+PG.vertex(2,1))/2-dy/10*PG.normal(2,end),['$e_{' num2str(k+1) '}$'],'interpreter','latex','fontsize',12)

        % vertices
%         plot(PG.vertex(1,:),PG.vertex(2,:),'ob','markerFaceColor','k','markersize',5)
%         text(PG.vertex(1,1)-0.7*mean(PG.normal(1,[1,PG.nv]))-0.5,PG.vertex(2,1)-0.7*mean(PG.normal(2,[1,PG.nv])),['s=' num2str(round(PG.sv(1),2))],'interpreter','latex')        
%         text(PG.vertex(1,1)-0.7*mean(PG.normal(1,[1,PG.nv]))-1.5,PG.vertex(2,1)-0.7*mean(PG.normal(2,[1,PG.nv]))-0.5,['(s=' num2str(round(PG.sv(end),2)) ')'],'interpreter','latex')        
% 
%         for k=2:PG.nv
%             text(PG.vertex(1,k)-0.7*mean(PG.normal(1,k-1:k))-0.5,PG.vertex(2,k)-0.7*mean(PG.normal(2,k-1:k)),['s=' num2str(round(PG.sv(k),2))],'interpreter','latex')
%         end
        text(0,-2,'$O_i$','interpreter','latex','fontSize',20)        
        text(0,-2,'$g$','interpreter','latex','fontSize',16,'color','r')  
        if ~isempty(PG.com)
                plot(PG.com(1,1),PG.com(2,1),'.r','markerFaceColor','r','markersize',20)
                quiver(PG.com(1,1),PG.com(2,1),PG.gv(1),0.5*PG.gv(2),'r','lineWidth',2)
%                 text(PG.com(1,1)+0.2,PG.com(2,1)-0.2,'CoM')
        end
%         grid on
        
        xl = [min(PG.vertex(1,:))-0.5*dx max(PG.vertex(1,:))+0.5*dx];
        yl = [min(PG.vertex(2,:))-0.5*dy max(PG.vertex(2,:))+0.5*dy];
        
%         text(0.9*xl(1)-0.1*dx,0.9*yl(2)-dy/10,'g','interpreter','latex','fontsize',10)
        hold off
  
        axis([xl yl]);
        axis equal
    end