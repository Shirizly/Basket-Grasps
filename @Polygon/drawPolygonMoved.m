function drawPolygonMoved(PG,axle,theta,d,color)
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        vert = PG.vertex;
        for i=1:size(PG.vertex,2)
            vert(:,i) = axle + R*(PG.vertex(:,i)-axle)+d;
        end
        hold on
        % edges
        plot(vert(1,:),vert(2,:),'-','Color',color)
        plot(vert(1,[1,end]),vert(2,[1,end]),'-','Color',color)
        
        % vertices
%         plot(vert(1,:),vert(2,:),'ob','markerFaceColor','k','markersize',5)
        
        
        if ~isempty(PG.com)
            comloc = R*(PG.com-axle)+axle+d;
                plot(comloc(1,1),comloc(2,1),'.r','markerFaceColor','r','markersize',15)
%                 quiver(comloc(1,1),comloc(2,1),0,-2,'r','maxHeadSize',10)
%                 text(comloc(1,1)-1.3,comloc(2,1)-0.5,'g','fontname','timesnewroman','color','r','fontsize',18)
        end 
%         grid on
%         xl = [min(vert(1,:))-2 max(vert(1,:))+2];
%         yl = [min(vert(2,:))-2 max(vert(2,:))+2];
%         
%         hold off
% 
%         axis([xl yl]);
        axis equal
    end
    