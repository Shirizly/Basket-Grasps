function drawPolygon2fGraspTable(PG,f1,f2,theta)
% s1 = 1.184;
% s2 = [];
% theta = -1.809;% -pi;
% f2 = [0;0];
% f1 = [-0.2;-0.03];
figure

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];


PG.drawPolygon();
plot(f1(1),f1(2),'.k','markerSize',12)
plot(f2(1),f2(2),'.k','markerSize',12)
hold on
fingerline = [f1,f2];
middlepoint = mean(fingerline,2);
direct = (fingerline-repmat(middlepoint,[1,2]))/0.9;
quiver(middlepoint(1)*[1,1],middlepoint(2)*[1,1],direct(1,:),direct(2,:),'k')
text(middlepoint(1),middlepoint(2)+0.2,'$\sigma$','Interpreter','latex','fontsize',18)
plot([middlepoint(1) middlepoint(1)+1],[middlepoint(2) middlepoint(2)],'k--')
text(middlepoint(1)+0.2,middlepoint(2),'$\phi$','Interpreter','latex','fontsize',18)
end