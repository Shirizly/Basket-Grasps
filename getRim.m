function [cont] = getRim(xx,tt,yy,U)
        % receives a specific inter finger distance (sig)
        % the inter-finger distance matrix Sigma
        % Outputs: cont - cell array of the contour lines of
        % the grasp, separated into the discontinuous parts

        figure
        [c,~] = contour(xx,tt,yy,[U,U],'showtext','on','linewidth',2,'linecolor','k');
        close
        % remove the out of bounds values
        c(:,(c(2,:)>max(tt))) = [];
        % separate into segments
        difc(1,:) = diff(c(1,:));
        difc(2,:) = diff(c(2,:));
        normdifc = sqrt(difc(1,:).*difc(1,:)+difc(2,:).*difc(2,:));
        sep = find(normdifc>0.1);
        nparts = length(sep)+1;
        sep = [0,sep,size(c,2)];
        part = cell(nparts,1);
        epsi = 0.5;
        mindist = 100;
        closest = nan;
        for i =1:nparts
            part{i} = c(:,sep(i)+1:sep(i+1));
            dist = min(diag(part{i}.'*part{i}));
            if dist<mindist
                mindist = dist;
                closest = i;
            end
        end
        % Identifying the particular contour for the grasp - tbd
         cont = part{closest};
    end