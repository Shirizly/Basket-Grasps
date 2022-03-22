function fig = plotGraphNodes(PG,cont,finger,nodes,otherFingerRelative,ContactHeight)
%% plot the COM height level sets
theta = -pi():0.1:pi();
Srange = PG.S;%PG.S(PG.S>=7&PG.S<=9.5);
[ height_matrix ] = single_support_height(PG,Srange,theta);
height_matrix = height_matrix + ContactHeight;
fig = figure;
ax(1)=axes;
hold on
c_layers=50;
 contour(Srange,theta,height_matrix.',c_layers,'linewidth',3)
 set(ax(1),'FontSize',18,'FontName','TimesNewRoman')
 xticks([0:2:PG.S(end)])
 xlabel(['$s_' num2str(finger) '$'],'fontSize',42,'interpreter','latex')
 ylabel('$\theta$','fontSize',42,'interpreter','latex')
h=colorbar;
set(get(h,'label'),'string','$\:U$','interpreter','latex','FontSize',40,'rotation',0);

%% plot ds contours  

for i =1:numel(cont)
    stmp = cont{i}(finger,:);
    thetatmp = cont{i}(4,:);
    plot(stmp,thetatmp,'k.','lineWidth',1)
end
%% fill DS penetration zones:
fill = 0;
if fill
    thetav = -pi()/2:0.01:pi()/2;
    
    Distance = findInnerPoints(PG,otherFingerRelative,thetav);
    ax(2) = axes;
    hold on
    contourf(PG.S,thetav,-100*Distance.',[0,0],'k')
    
    contour(PG.S,thetav,-100*Distance.',[0,0],'k')
    colormap(ax(2),[0 0 0])
    %move x- and y-axis of axes2
    pos=get(ax(1),'Position');
    set(ax,'Position',pos);
    set(ax(2), 'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XTick',[],...
        'YTick',[]);
end
%% mark edges of strips (and half strips) in dashes
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
MS = 20;%10;
LW = 6;%4;
if ~isempty(ss_hook)
hook = ss_hook(ss_hook(:,5)==finger,:); %hooks
plot(hook(:,1),hook(:,3),'o','color',[0,0.9,0],'markerSize',MS,'lineWidth',LW)
end
saddle = ss_saddle(ss_saddle(:,4)==finger,:); % mark ss saddle point equilibrium
plot(saddle(:,1),saddle(:,3),'ob','markerSize',MS,'lineWidth',LW)

ssmin = ss_min(ss_min(:,5)==finger,:); % mark ss vertex minimum points
plot(ssmin(:,1),ssmin(:,3),'o','color',[0,0.9,0],'markerSize',MS,'lineWidth',LW)

ssmax = ss_max(ss_max(:,4)==finger,:); % mark ss vertex minimum points
plot(ssmax(:,1),ssmax(:,3),'o','color',[1,0,0],'markerSize',MS,'lineWidth',LW)
%% double support markers
virtual = ds_virtual; % mark ds virtual nodes (on edges of strips)
plot(virtual(:,finger),virtual(:,4),'+b','markerSize',8,'lineWidth',3)

dsmax = ds_max; % mark ds max nodes
plot(dsmax(:,finger),dsmax(:,4),'+r','markerSize',8,'lineWidth',3)

dsmin = ds_min; % mark ds min nodes
plot(dsmin(:,finger),dsmin(:,4),'+','color',[0,0.9,0],'markerSize',8,'lineWidth',3)

%% remove markers protruding out of the axes
% xl = get(gca,'XLim');
% yl = get(gca,'YLim');
% set(gca,'clipping','off')
% extremes = [xl(2)-xl(1), yl(2)-yl(1)];
% rectangle('Position',[xl(1)-extremes(1), yl(2)            , 3*extremes(1),   extremes(2)],'FaceColor',[1 1 1],'EdgeColor','none'); % Top
% rectangle('Position',[xl(1)-extremes(1), yl(1)-extremes(2), 3*extremes(1),   extremes(2)],'FaceColor',[1 1 1],'EdgeColor','none'); % Bottom
% rectangle('Position',[xl(2)            , yl(1)-extremes(2),   extremes(1), 3*extremes(2)],'FaceColor',[1 1 1],'EdgeColor','none'); % Right
% rectangle('Position',[xl(1)-extremes(1), yl(1)-extremes(2),   extremes(1), 3*extremes(2)],'FaceColor',[1 1 1],'EdgeColor','none'); % Left
% set(gca,'XLim',xl);
% set(gca,'YLim',yl);
% set(gca,'box','on')
% set(gca,'Layer','top')

hold off
end

