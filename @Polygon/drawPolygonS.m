function fig = drawPolygonS(PG,s1,s2,linesegs)
        % draw the poygon and add theline segments representing possible
        % grasps
        fig = PG.drawPolygon();
        hold on
        for i=1:3:length(s1)
            sind1 = find(PG.S==s1(i));
            sind2 = find(PG.S==s2(i));
            plot(PG.X(sind1,1),PG.X(sind1,2),'.c','markerSize',6);
            plot(PG.X(sind2,1),PG.X(sind2,2),'.c','markerSize',6); 
        end
        for i=1:size(linesegs,1)
            sind1 = find(PG.S==linesegs(i,1));
            sind2 = find(PG.S==linesegs(i,2));
            sind3 = find(PG.S==linesegs(i,3));
            sind4 = find(PG.S==linesegs(i,4));
            plot(PG.X([sind1,sind3],1),PG.X([sind1,sind3],2),'.c','markerSize',6);
            plot(PG.X([sind2,sind4],1),PG.X([sind2,sind4],2),'c','linewidth',3); 
        end
        hold off 
    end
    