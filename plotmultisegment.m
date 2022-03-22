function plotmultisegment(varargin)
rate = [0.5,0.27,0.2];
X = []; Y1 = []; Y2 = []; Y3 = [];
for i=1:nargin
    filename = varargin{i};
    load(filename,'BGDepth2','BGDepth','SSDepth');
    BGS = linspace(0,1,size(SSDepth,2));
%     BGSn = linspace(0,0.98,size(BGdepthlist,2));
    % figure;
    hold on;
    correc = sum(rate(1:i-1))*1.05;
    BGu = BGS*rate(i);

        X = [X,correc+BGu];
        Y1 = [Y1,BGDepth2];
        Y2 = [Y2,BGDepth];
        Y3 = [Y3,SSDepth];
    
end
p1 = plot(X,Y1,'g','lineWidth',3,'DisplayName','depth graph');
% p2 = plot(correc+BGSn*rate(i),BGdepthlist,'+k','MarkerSize',8,'lineWidth',3,'DisplayName','graph search');
p3 = plot(X,Y2,'r--','lineWidth',2,'DisplayName','double-support');
p4 = plot(X,Y3,'b--','lineWidth',2,'DisplayName','single-support');
% title('Basket Depth for Basket grasp segment - Cup','interpreter','latex','fontSize',14);
xlabel('basket grasp curve parameter $u$','interpreter','latex');
ylabel('$\Delta U$','interpreter','latex');
legend([p1,p3,p4])
hold off
end