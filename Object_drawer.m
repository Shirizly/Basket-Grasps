clear all
close all

%%
% load object.mat
% object = [object(:,7:8) object(:,1:6)];
% object = [object(:,1:3) [-0.15;0.23] [-0.1;0.21] object(:,4:8)];
% PG = Polygon(object,0*com);
% polyDrawing = PG.drawPolygon();

i=1;
texthandle = [];
flag = true;
figure
grid on
hold on
% imshow([pwd '\Examples for research proposal\Gun\snip.png'])
hold on
set(gca,'visible',1)
% xlim = get(gca,'xlim');
% ylim = get(gca,'ylim');
axis([-5 5 -5 5])
xv = [min(xlim),min(xlim)+diff(xlim)/10,min(xlim)+diff(xlim)/10,min(xlim)];
yv = [min(ylim),min(ylim),min(ylim)+diff(ylim)/10,min(ylim)+diff(ylim)/10];
pos = [xv(1),yv(1),max(xv)-min(xv),max(yv)-min(yv)];
rectangle('Position',pos,'faceColor',[1 1 0 0.2])
texthandle = text(mean(xv),mean(yv),'complete')


while flag
%     if ~isempty(texthandle)
%         delete(texthandle);
%     end
%     texthandle = text(0,5,num2str(i));
    [object_x(i),object_y(i)] = ginput(1);
    
    if inpolygon(object_x(i),object_y(i),xv,yv)
        flag = false;
        object_x(end) = [];
        object_y(end) = [];
        break;
    end
    plot(object_x(i),object_y(i),'.k');
    if i>1
      plot(object_x(i-1:i),object_y(i-1:i),'k'); 
    end

    i = i+1;
end
plot([object_x(1),object_x(end)],[object_y(1),object_y(end)],'k'); 
delete(texthandle);
texthandle = text(0,5,'center of mass');
[com(1),com(2)] = ginput(1);
% plot(com(1),com(2),'*r');

object = [object_x;object_y]
com = com.'

%%
Round_object = round(object,2);
% com(2) = max(Round_object(2,:))-com(2);
% Round_object(2,:) = (max(Round_object(2,:))-Round_object(2,:));
figure
hold on
plot(Round_object(1,:),Round_object(2,:))
plot(com(1),com(2),'+')
axis equal
%%

save('objectHiRes1.mat','Round_object','com')
