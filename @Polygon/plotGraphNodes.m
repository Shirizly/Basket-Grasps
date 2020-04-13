function fig = plotGraphNodes(PG,cont,finger,nodes,otherFingerRelative,ContactHeight)
%% plot the COM height level sets
theta = -pi():0.01:pi();
[ height_matrix ] = single_support_height(PG,PG.S,theta);
height_matrix = height_matrix + ContactHeight;
fig = figure;
hold on
c_layers=50;
contour(PG.S,theta,height_matrix.',c_layers)
xlabel(['s_' num2str(finger)],'fontSize',22)
ylabel('\theta','fontSize',22)
colorbar
%% plot ds contours and strip edges 
for i =1:numel(cont)
    stmp = cont{i}(finger,:);
    thetatmp = cont{i}(4,:);
    plot(stmp,thetatmp,'k.','lineWidth',3)
end

%% fill DS penetration zones: (only to make pretty figures, takes time and ladens computer)
% thetav = -pi():0.02:pi();
% Distance = findInnerPoints(PG,otherFingerRelative,thetav);
% contour(PG.S,thetav,Distance.',[-10:0.01:0],'showtext','off','linewidth',4,'linecolor','k')
%% mark edges of strips (and half strips) in 
for i= 2:length(PG.vertex(1,:)) % mark edges and feasibility edges of theta
    plot(PG.S(PG.VL(i))*ones(2,1),[-pi(),pi()],'--k')
    theta0(1) = wrapToPi(pi()/2-atan2(PG.normal(2,i),PG.normal(1,i)));
    theta0(2) = wrapToPi(theta0(1)+pi());
    plot(PG.S(PG.VL([i,i+1])),theta0(1)*ones(2,1),'--k');
    plot(PG.S(PG.VL([i,i+1])),theta0(2)*ones(2,1),'--k');

%     %% adding rectangular markers for unfeasible support areas
%     theta0(1) = wrapToPi(pi()/2-atan2(PG.normal(2,i),PG.normal(1,i)));
%     theta0(2) = wrapToPi(theta0(1)+pi());
%     if theta0(1)>=0
%         pos = [PG.S(PG.VL(i)),theta0(2),PG.S(PG.VL(i+1))-PG.S(PG.VL(i)),theta0(1)-theta0(2)];
%         rectangle('Position',pos,'faceColor',[0 0 0 0.1])
%     else
%         pos1 = [PG.S(PG.VL(i)),theta0(2),PG.S(PG.VL(i+1))-PG.S(PG.VL(i)),pi()-theta0(2)];
%         pos2 = [PG.S(PG.VL(i)),-pi(),PG.S(PG.VL(i+1))-PG.S(PG.VL(i)),pi()+theta0(1)];
%         rectangle('Position',pos1,'faceColor',[0 0 0 0.1])
%         rectangle('Position',pos2,'faceColor',[0 0 0 0.1])
%     end
end


ds_max =  nodes{1};
ds_min= nodes{2};
ds_virtual = nodes{3};
ss_max = nodes{4};
ss_min= nodes{5};
ss_saddle = nodes{6};
ss_hook = ss_min(find(ss_min(:,4)==1),:); %#ok<*FNDSB>
%% single support markers
hook = ss_hook; %hooks
plot(hook(:,1),hook(:,3),'o','color',[0,0.9,0],'markerSize',10,'lineWidth',4)

saddle = ss_saddle(ss_saddle(:,4)==finger,:); % mark ss saddle point equilibrium
plot(saddle(:,1),saddle(:,3),'ob','markerSize',10,'lineWidth',2)

ssmin = ss_min(ss_min(:,5)==finger,:); % mark ss vertex minimum points
plot(ssmin(:,1),ssmin(:,3),'o','color',[0,0.9,0],'markerSize',8,'lineWidth',2)

ssmax = ss_max(ss_max(:,4)==finger,:); % mark ss vertex minimum points
plot(ssmax(:,1),ssmax(:,3),'o','color',[1,0,0],'markerSize',8,'lineWidth',2)
%% double support markers
virtual = ds_virtual; % mark ds virtual nodes (on edges of strips)
plot(virtual(:,finger),virtual(:,4),'+b','markerSize',8,'lineWidth',3)

dsmax = ds_max; % mark ds max nodes
plot(dsmax(:,finger),dsmax(:,4),'+r','markerSize',8,'lineWidth',3)

dsmin = ds_min; % mark ds min nodes
plot(dsmin(:,finger),dsmin(:,4),'+','color',[0,0.9,0],'markerSize',8,'lineWidth',3)

hold off
end