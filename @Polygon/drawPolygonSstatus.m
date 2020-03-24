function fig = drawPolygonSstatus(PG,s1,s2,linesegs,status)
        % draw the poygon and add the each possible grasp
        fig = PG.drawPolygon();
        hold on
        col = [0 1 0; 0 0 0; 1 0 0];
        for i=1:size(linesegs,1)
            sind1 = find(PG.S==linesegs(i,1));
            sind2 = find(PG.S==linesegs(i,2));
            sind3 = find(PG.S==linesegs(i,3));
            sind4 = find(PG.S==linesegs(i,4));
            plot(PG.X([sind1,sind3],1),PG.X([sind1,sind3],2),'.k','markerSize',4);
            plot(PG.X([sind2,sind4],1),PG.X([sind2,sind4],2),'k','linewidth',2); 
        end
        
        for i=1:3:length(s1)
            sind1 = find(PG.S==s1(i));
            sind2 = find(PG.S==s2(i));
            plot(PG.X(sind1,1),PG.X(sind1,2),'.','markerSize',10,'MarkerEdgeColor',col(status(i)+2,:));
            plot(PG.X(sind2,1),PG.X(sind2,2),'.','markerSize',10,'MarkerEdgeColor',col(status(i)+2,:)); 
        end
        hold off 
    end
    