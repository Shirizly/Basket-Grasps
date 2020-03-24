function fig = drawPolyGrasp(PG,sf1,sf2,f1pos,f2pos)
fig = PG.drawPolygon();
hold on
R = 0.4;
n1 = R*PG.get('normal',sf1);
n2 = R*PG.get('normal',sf2);
%         t = linspace(0,2*pi);
cent1 = f1pos-n1;
plot(f1pos(1),f1pos(2),'.k','markersize',10)
text(cent1(1)-0.5*R,cent1(2),'$O_1$','fontsize',12,'interpreter','latex')
% plot(cent1(1)+R*cos(t),cent1(2)+R*sin(t),'k')
cent2 = f2pos-n2;
if ~isempty(f2pos)
    plot(f2pos(1),f2pos(2),'.k','markersize',10)
    text(cent2(1)-0.5*R,cent2(2),'$O_2$','fontsize',12,'interpreter','latex')
end
% plot(cent2(1)+R*cos(t),cent2(2)+R*sin(t),'k')
end
    