clc
clearvars 
close all
contspace = [];
distmag = 0.2;
R = 5;
beta = 2*pi()/10;



pwd = 'Examples for research proposal\double support\random4\';
nsample = 50;
dimv = 16:4:64;
maxv = zeros(length(dimv),1);
meanmaxv = zeros(length(dimv),1);
meanv = zeros(length(dimv),1);
%%
% for j=8:length(dimv)

for j = 5
    dim = dimv(j)
    object = zeros(2,dim);
    tet = linspace(0,2*pi(),dim+1);
    tet(end) = [];
    for l = 1:nsample
        for k=1:dim
            object(:,k) = (0.1+rand(1))*R*[cos(tet(k));sin(tet(k))];
        end
        com = [0;0];
        PG = Polygon(object,com);
        PG.drawPolygonClear();
        title(['n = ' num2str(dim)],'interpreter','latex')
        saveas(gcf,[pwd 'comb' num2str(j) ',' num2str(l) '.bmp'])
        close(gcf)
        res = dim;
        
        [PG,S,X,VL] = PG.findBdyVariable(res);
        contspace{end+1} = figure;
        hold on
        for i=1:PG.nv
            plot([0 PG.S(end)],PG.S(PG.VL(i))*ones(1,2),'k--','linewidth',1)
            plot(PG.S(PG.VL(i))*ones(1,2),[0 PG.S(end)],'k--','linewidth',1)
        end
        xlim([0 PG.S(end)]);
        ylim([0 PG.S(end)]);
%         c_layers=10;
%         Sigma=inter_finger_distance(X,X);
%         contour(S,S,Sigma,c_layers)
        xlabel('$s_1$','interpreter','latex','fontsize',18)
        ylabel('$s_2$','interpreter','latex','fontsize',18)
%         colorbar
        tic
        N=40;
        numsegvector = zeros(N,1);
        count = 1;
        for i=1:N/2
            s1init = (rand(1))/2*PG.S(end);
            s2init = (1+rand(1))/2*(PG.S(end));
%             s1init = i/N*PG.S(end)/4+rand(1)*1E-5;
%             s2init = i/N*PG.S(end)-rand(1)*1E-5;
%             s1init = rand(1)*2E-1;
%             s2init = PG.S(PG.VL(k))-rand(1)*(PG.S(PG.VL(k))-PG.S(PG.VL(k-1)));
            plot(s1init,s2init,'+r')
            numsegvector(count) = equisigmacontourNum(PG,s1init,s2init);
            if numsegvector(count) > 20*dim
                susp = numsegvector(count)
                i = i-1;
                count = count-1;
            end
            count = count+1;
        end
        for i=N/2+1:N
            itag = ceil((i-N/2)/(N/2)*PG.nv)+1;
            s1init = rand(1)*2E-1;
            s2init = PG.S(PG.VL(itag))-rand(1)*(PG.S(PG.VL(itag))-PG.S(PG.VL(itag-1)));
            plot(s1init,s2init,'+r')
            numsegvector(count) = equisigmacontourNum(PG,s1init,s2init);
            if numsegvector(count) > 20*dim
                susp = numsegvector(count)
                i = i-1;
                count = count-1;
            end
            count = count+1;
        end
        time = toc
        maxv(j) = maxv(j)+ max(numsegvector);
        title(['n = ' num2str(dim)],'interpreter','latex')
        saveas(gcf,[pwd 'comb' num2str(j) ',' num2str(l) ' curves.bmp'])
        close(gcf)
        meanv(j) = mean(numsegvector)
    end
    meanmaxv(j) = maxv(j)/nsample
end
%%
figure
plot(dimv,meanmaxv)
title('Randomly sampled equi-\sigma curves')
xlabel('n - number of object vertices')
ylabel('maximum segment count')
savefig(gcf,[pwd ,'Randomly sampled curves.fig'])
saveas(gcf,[pwd ,'Randomly sampled curves.bmp'])
%%
figure
plot(dimv,meanmaxv.'./dimv)
title('Randomly sampled equi-\sigma curves - normalized')
xlabel('n - number of object vertices')
ylim([0,10])
ylabel('maximum segment count, normalized')
savefig(gcf,[pwd ,'Randomly sampled curves normalized.fig'])
saveas(gcf,[pwd ,'Randomly sampled curves normalized.bmp'])
%% Fit a linear/quadratic function
xy = [dimv.',meanmaxv];
b = polyfit(xy(:,1), xy(:,2), 2)                  % Estimate Parameters For Quadratic Fit
y_fit = polyval(b, xy(:,1));                        % Evaluate Fitted Curve
figure(1)
plot(xy(:,1), xy(:,2), 'pg')
hold on
plot(xy(:,1), y_fit, '-r')
hold off
grid
sum(((y_fit-xy(:,2))./xy(:,2)).^2)
%% object constructors:
% square double-sided comb:
% corner1 = [-5;-5];
% corner2 = [5;5];
% corner3 = [-5;5];
% corner4 = [5;5];

% object(1,1:dim/2) = linspace(corner1(1),corner2(1),dim/2);
% object(1,dim/2+1:end) = linspace(corner2(1),corner1(1),dim/2);
% object(2,1:2:dim/2) = corner1(2)*ones(1,dim/4);
% object(2,2:2:dim/2) = (1-distmag)*corner1(2)*ones(1,dim/4);
% object(2,dim/2+1:2:end) = corner3(2)*ones(1,dim/4);
% object(2,dim/2+2:2:end) = (1+distmag)*corner3(2)*ones(1,dim/4);

% circular comb:
% distmag = 0.2;
% R = 5;
% beta = 2*pi()/10;

% object = zeros(2,dim);
% tet = linspace(0,2*pi(),dim+1);
% tet(end) = [];
% object(:,1:2:end) = R*[cos(tet(1:2:end));sin(tet(1:2:end))];
% fix = R*tan(2*pi()/dim)/tan(beta);
% object(:,2:2:end) = (1+fix)*R*[cos(tet(2:2:end));sin(tet(2:2:end))];