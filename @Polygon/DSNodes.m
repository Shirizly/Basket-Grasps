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
%         difh = diff(dsh);
%         if difh(1)<0
%             maxinds_temp = [1,maxinds_temp];
%         else
%             mininds_temp = [1,mininds_temp];
%         end
%         if difh(end)>0
%             maxinds_temp = [maxinds_temp,length(dsh)];
%         else
%             mininds_temp = [mininds_temp,length(dsh)];
%         end
        for i = 1:length(maxinds_temp)
            ds_max(size(ds_max,1)+1,:) = [dss1(maxinds_temp(i)),dss2(maxinds_temp(i)),dsh(maxinds_temp(i)),dstheta(maxinds_temp(i)),j,maxinds_temp(i)];
        end
        for i = 1:length(mininds_temp)
            ds_min(size(ds_min,1)+1,:) = [dss1(mininds_temp(i)),dss2(mininds_temp(i)),dsh(mininds_temp(i)),dstheta(mininds_temp(i)),j,mininds_temp(i)];
        end
    end
    
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