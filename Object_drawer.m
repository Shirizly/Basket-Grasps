clear all
close all

%%
i=1;
texthandle = [];
flag = true;
figure
grid on
hold on
xv = [-5,-3,-3,-5];
yv = [0,0,1,1];
pos = [xv(1),yv(1),max(xv)-min(xv),max(yv)-min(yv)];
rectangle('Position',pos,'faceColor',[0 0 0 0.2])
text(mean(xv),mean(yv),'complete')
axis([-5,5,-5,5])
while flag
    if ~isempty(texthandle)
        delete(texthandle);
    end
    texthandle = text(0,5,num2str(i));
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
plot(com(1),com(2),'*r');

object = [object_x;object_y]
com = com.'

%%
Round_object = round(object,1)
% figure
% plot(Round_object(1,:),Round_object(2,:))
%%
save('object.mat','object','com')