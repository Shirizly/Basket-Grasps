function [ds_max,ds_min,ds_virtual] = DSNodes(PG,cont)
%function receives Polygon object and full double-support polygon
% (cell array of segments of continuous double-support polygon [s1;s2;h;theta;d]
%function outputs lists of double-support nodes for the graph, including
%indexs in the contour,and parameters s,h,theta.

ds_max = [];
ds_min = [];
ds_virtual = [];
nseg = numel(cont);

epsi = max(abs(diff(PG.S)));

for j = 1:nseg
    contv = cont{j};
    dss1 = contv(1,:);
    dss2 = contv(2,:);
    dstheta = contv(4,:);
    dsh = contv(3,:);
    if length(dsh)>2
        [~,maxinds_temp] = findpeaks(dsh);
        [~,mininds_temp] = findpeaks(-dsh);
        
        % remove nodes that are on top of each other (
        removemax = []; removemin = [];
        for i = 1:length(maxinds_temp) 
            identical = find(abs(mininds_temp-maxinds_temp(i))<=1);
            if ~isempty(identical)
                removemax = [removemax i];
                removemin = [removemin identical];
            end
        end
        maxinds_temp(removemax) = [];
        mininds_temp(removemin) = [];
%         remove = [];
%         for i = 2:length(maxinds_temp) 
%             between = find(mininds_temp<maxinds_temp(i) & mininds_temp>maxinds_temp(i-1));
%             if isempty(between)
%                 remove = [remove i]
%             end
%         end
%         maxinds_temp(remove) = [];
%         remove = [];
%         for i = 2:length(mininds_temp) 
%             between = find(maxinds_temp<mininds_temp(i) & maxinds_temp>mininds_temp(i-1)));
%             if isempty(between)
%                 remove = [remove i]
%             end
%         end
%         mininds_temp(remove) = [];

        for i = 1:length(maxinds_temp)
            ds_max(size(ds_max,1)+1,:) = [dss1(maxinds_temp(i)),dss2(maxinds_temp(i)),dsh(maxinds_temp(i)),dstheta(maxinds_temp(i)),j,maxinds_temp(i)];
        end
        for i = 1:length(mininds_temp)
            ds_min(size(ds_min,1)+1,:) = [dss1(mininds_temp(i)),dss2(mininds_temp(i)),dsh(mininds_temp(i)),dstheta(mininds_temp(i)),j,mininds_temp(i)];
        end
    end
    %% marking the virtual nodes that are max/min properly (started working on this then realized it doesn't matter now)
%     end_connector = [];
%     range = [1:nseg];
%     range(j) = [];
%     % check if the segment is circular:
%     if dss1(end)==dss1(1) && dss2(end) == dss2(1)
%         end_connector = [j,1];
%     end
%     
%     
%     for k = range
%         if ~isempty(end_connector)
%             break;
%         end
%         contv2 = cont{k};
%         dss12 = contv(1,:);
%         dss22 = contv(2,:);
%         if dss1(end)==dss12(1) && dss2(end) == dss22(1)
%             end_connector = [k,1];
%         end
%         if dss1(end)==dss12(end) && dss2(end) == dss22(end)
%             end_connector = [k,length(dss22)];
%         end
%     end
%     % having found the end:
%     contv2 = cont{end_connector(1)};
%         dss12 = contv(1,:);
%         dss22 = contv(2,:);
%         dsh2 = contv(3,:);
%         
    
    %% adding virtual nodes every time the contour crosses a vertex with either finger
    for k=1:length(PG.VL)
        diff_VL = abs(dss1-PG.S(PG.VL(k)));
        [~,find_start] = findpeaks(-diff_VL);
        find_start = [1,find_start,length(dss1)];
        for i = 1:length(find_start)
            if diff_VL(find_start(i))<epsi
                ds_virtual(size(ds_virtual,1)+1,:) = [dss1(find_start(i)),dss2(find_start(i)),dsh(find_start(i)),dstheta(find_start(i)),j,find_start(i)];
            end
        end
        
        diff_VL = abs(dss2-PG.S(PG.VL(k)));
        [~,find_start] = findpeaks(-diff_VL);
        find_start = [1,find_start,length(dss2)];
        for i = 1:length(find_start)
            if diff_VL(find_start(i))<epsi
                ds_virtual(size(ds_virtual,1)+1,:) = [dss1(find_start(i)),dss2(find_start(i)),dsh(find_start(i)),dstheta(find_start(i)),j,find_start(i)];
            end
        end
        
        %% calculating the saddle-angles, adding virtual nodes whenever the contour crosses them:
        if k<length(PG.VL)
            normal = PG.normal(:,k);
            saddleTheta1 = wrapToPi(pi()/2-atan2(normal(2),normal(1))); %angle of the higher saddle line
            diff_theta1 = abs(dstheta-saddleTheta1);
            sStart = PG.S(PG.VL(k));
            sEnd = PG.S(PG.VL(k+1));
            nodesINstrip = (sStart<dss1&dss1<sEnd)|(sStart<dss2&dss2<sEnd);
            [~,find_cross] = findpeaks(-diff_theta1);
            find_cross(find(nodesINstrip(find_cross)==0)) = [];
            for i = 1:length(find_cross)
                if diff_theta1(find_cross(i))<epsi
                    ds_virtual(size(ds_virtual,1)+1,:) = [dss1(find_cross(i)),dss2(find_cross(i)),dsh(find_cross(i)),dstheta(find_cross(i)),j,find_cross(i)]; %#ok<*AGROW>
                else
%                     disp(diff_theta1(find_cross(i)))
                end
            end            
            saddleTheta2 = wrapToPi(pi()+saddleTheta1); %angle of the higher saddle line
            diff_theta2 = abs(dstheta-saddleTheta2);
            [~,find_cross] = findpeaks(-diff_theta2);
            find_cross(find(nodesINstrip(find_cross)==0)) = []; %#ok<*FNDSB>

            for i = 1:length(find_cross)
                if diff_theta2(find_cross(i))<epsi
                    ds_virtual(size(ds_virtual,1)+1,:) = [dss1(find_cross(i)),dss2(find_cross(i)),dsh(find_cross(i)),dstheta(find_cross(i)),j,find_cross(i)]; %#ok<*AGROW>
                else
%                     disp(diff_theta2(find_cross(i)))
                end
            end 
            
            

        end
    end 
    
end
ds_virtual = remove_minmax(ds_virtual,[ds_min;ds_max]);
ds_virtual = remove_copies(ds_virtual);
end

function list = remove_copies(list)
remove = [];
for i=1:size(list,1)
    for j=i+1:size(list,1)
        if list(i,:)==list(j,:)
            remove = [remove;j];
        end
    end
end
list(remove,:) = [];
end

function list = remove_minmax(list,list2)
remove = [];
for i=1:size(list,1)
    for j=1:size(list2,1)
        if list(i,:)==list2(j,:)
            remove = [remove;i];
        end
    end
end
list(remove,:) = [];
end