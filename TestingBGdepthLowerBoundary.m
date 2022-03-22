clearvars
%% Define the polygon
% %% Define the polygon
clc

object = 11;
symmetry = 0;
switch object
    case 0 % pentagon
        theta = 0;%-pi()/2;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
        PG = Polygon(R*[3,-4;3,4;-1.5,5;-3,0;-1.5,-5].'+[0;3],[0;3])%0.75,-4.5
        subFolderName = 'Examples for research proposal\Pentagon'; % where to save results
        res = 50;
        symmetry = [1,2,3,3,2;1,5,4,3,2;-1,-1,-5,-4,-5;-2,-1,-5,-4,-3];
    case 1 % cup
       PG = Polygon([-2,-2;2,-2;2,1;4,1;4,2;-3,2].'+[0;3],[0;3]);
       subFolderName = 'Examples for research proposal\Cup'; % where to save results

%        basepos = [f1,f2];
       res = 50;
    case 2 % bottle
        load object5.mat
        %     object = [object(:,7:8) object(:,1:6)];
%         com = [0.045,0.07].';
%         object = object + 0.001*randn(size(object));
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Bottle'; % where to save results
        f2 = [0.06;0.01];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 1500;
        starting_node = [2,21];
    case 3 % Elon
        load object.mat
        %     object = [object(:,7:8) object(:,1:6)];
        object = Round_object;
        com = [0.045,0.07].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Elon'; % where to save results
    case 4 % Crown
        scale = 0.02;
        scale = 1;
        object = scale*[-2,-1;2,-1;3,3;1,0;0,1;-1,0;-4,1].';
        com = [0,-0.00].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Crown'; % where to save results
        f2rel = [2.7519;2.0074];
        f1rel = [-2.8470;-0.1530];
        f1 = [0;0];
        f2 = f2rel-f1rel;
        res = 50;
         starting_node = [2,4];    
    case 5 % new bottle
        load NewBottle.mat
        %     object = [object(:,7:8) object(:,1:6)];
%         com = [0.045,0.07].';
%         object = object + 0.001*randn(size(object));
%         theta =  1.515;
%         R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
%         object = R*object;
        com(2) = max(object(2,:))-com(2);
        object(2,:) = max(object(2,:))-object(2,:);
        object(2,1) = 0;
        object = round(object./20,1);
        com = round(com./20,0);
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Bottle'; % where to save results
        f2 = [4;2];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 50;
        starting_node = [2,24];
        case 6 % Gun
%         load gun2.mat
%         object = Round_object;
%         object = 0.8*[1,10;26,10;21,1;29,9;31,2;28,0;32,0;32,12;35,16;32,18;33,16;32,14;31,14;31,16;1,16;0.5,13.5;0,13;0.5,12.5].';
        
        object = 0.8*[1,10;26,10;21,1;29,9;32,0;32,12;36,12;40,8;41,9;35,14;31,16;1,16].';%32,14;31,14
        com = 0.8*[27;10];
        theta =  0.423;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
         object = R*object;
         com = R*com;
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Gun'; % where to save results
        f2 = [18.2;16.09];
        f1 = [11.4;13.9];
        basepos = [f1,f2];
        res = 30;
        starting_node = [2,6];
        
        case 7 % Crown with more edges
        object = [-2,-1;2,-1;2.25,1;3.25,3;1,0;0,1;-1,0;-4,1;-2.75,0].';
        com = [0,0].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Crown'; % where to save results
        f2rel = [2.7519;2.0074];
        f1rel = [-2.8470;-0.1530];
        f1 = [0;0];
        f2 = f2rel-f1rel;
        res = 50;
         starting_node = [2,4];  
         
    case 8 % New object for ss escape
        load objectConnector2.mat
         object = Round_object(:,end:-1:1);
%          object = [object(:,4:end),object(:,1:3)];
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\SSEscape'; % where to save results
        res = 150;
    case 9 % object with caging and immobilization
        scale = 0.02;
        object = scale*[1,-2;2,1;2,2;3,3;-3,3;-1,1.5;-2,-2].';
        com = [0,0].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\MultiPeak'; % where to save results
        res = 1500;
    case 10 % object with many vertices
        load objectHiRes1.mat
         object = Round_object(:,end:-1:1);
%          object = [object(:,4:end),object(:,1:3)];
        PG = Polygon(object,com);
        res = 50;
    case 11 %object to play with escape features
        object = [2,-3;4,2;5,3;-5,3;-2,-3].';
        com = [0,0].';
        PG = Polygon(object,com);
        res = 50;
end

polyDrawing = PG.drawPolygon();
hold off
%% Find equilibrium grasps
tic
[PG,S,X,VL] = PG.findBdyVariable(res);

[s1,s2,~] = PG.Eqcheck();
EqCurves = PG.EqCurveFinder();
for i=1:numel(EqCurves)
    curve = EqCurves{i};
    s1 = [s1, linspace(curve(1,1),curve(1,2),res)];
    s2 = [s2, linspace(curve(2,1),curve(2,2),res)];
    s1 = [s1, linspace(curve(2,1),curve(2,2),res)];
    s2 = [s2, linspace(curve(1,1),curve(1,2),res)];
end

% Check grasps for minimum/maximum
[BG,NBG] = PG.StatusSeparate(s1,s2,X);
timeSynesthise = toc
%% Full contact space figure
switch 0
    case 0
        contspace = figure;
hold on
c_layers=10;
Sigma=inter_finger_distance(X,X);
contour(S,S,Sigma,c_layers)
xlabel('$s_1$','interpreter','latex','fontsize',28)
ylabel('$s_2$','interpreter','latex','fontsize',28)
h = colorbar;
set(get(h,'label'),'string','$\sigma$','interpreter','latex','FontSize',28,'rotation',0);
col = [0 0.7 0; 0 0 0; 1 0 0];
% plot(BG(1,:),BG(2,:),'.','MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(1,:),'markersize',8)
% plot(NBG(1,:),NBG(2,:),'.','MarkerEdgeColor',col(3,:),'MarkerFaceColor',col(3,:),'markersize',8)
% for i=1:numel(EqCurves)
%     curve = EqCurves{i};
%     plot(curve(1,:),curve(2,:),'k--')
%     plot(curve(2,:),curve(1,:),'k--')
% end
for i=1:length(BG)
    plot(BG(1,i),BG(2,i),'.','MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(1,:),'markersize',8)
end
for i=1:length(NBG)
    plot(NBG(1,i),NBG(2,i),'.','MarkerEdgeColor',col(3,:),'MarkerFaceColor',col(1,:),'markersize',8)
end
for i=1:PG.nv
    plot([0 PG.S(end)],PG.S(PG.VL(i))*ones(1,2),'k--','linewidth',1)
    plot(PG.S(PG.VL(i))*ones(1,2),[0 PG.S(end)],'k--','linewidth',1)
end

% xl = [PG.S(PG.VL(5)) PG.S(PG.VL(7))];
% yl = [PG.S(PG.VL(2)) PG.S(PG.VL(4))];
% set(gca,'xlim',xl)
% set(gca,'ylim',yl)

if ~exist('BGSegEnds','var')
% Identify segment of basket grasps by pointing at end points:
e1 = ginput(1);
e2 = ginput(1);
distfrome1 = (BG(1,:)-e1(1)).^2+(BG(2,:)-e1(2)).^2;
inde1 = find(distfrome1==min(distfrome1),1,'first');
distfrome2 = (BG(1,:)-e2(1)).^2+(BG(2,:)-e2(2)).^2;
inde2 = find(distfrome2==min(distfrome2),1,'first');
BGSegEnds = [BG(1,inde1),BG(1,inde2);BG(2,inde1),BG(2,inde2)];
end
plot(BGSegEnds(1,:),BGSegEnds(2,:),'+','MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(1,:),'markerSize',10,'lineWidth',4)
set(contspace,'renderer','painters')
hold off
    case 1
        openfig('CrownBasketCurveNew.fig');
        if ~exist('BGSegEnds','var')
            % Identify segment of basket grasps by pointing at end points:
            e1 = ginput(1);
            e2 = ginput(1);
            distfrome1 = (BG(1,:)-e1(1)).^2+(BG(2,:)-e1(2)).^2;
            inde1 = find(distfrome1==min(distfrome1),1,'first');
            distfrome2 = (BG(1,:)-e2(1)).^2+(BG(2,:)-e2(2)).^2;
            inde2 = find(distfrome2==min(distfrome2),1,'first');
            BGSegEnds = [BG(1,inde1),BG(1,inde2);BG(2,inde1),BG(2,inde2)];
        end
end



%%
%%
drawPolygon2fGraspRange(PG,BGSegEnds);

%% find each depth curve iteratively, starting from the basket grasp feature pair:
tic

reso = 100;
BGS = linspace(0,1,reso);
BGP1 = zeros(2,reso);
BGP2 = zeros(2,reso);
BGPrel = zeros(2,reso);
BGPhi = zeros(1,reso);
BGs1 = linspace(BGSegEnds(1,1),BGSegEnds(1,2),length(BGS));
BGs2 = linspace(BGSegEnds(2,1),BGSegEnds(2,2),length(BGS));
% BGs1([1,end]) = [];
% BGs2([1,end]) = [];
BGSeg = [BGs1;BGs2];
for i=1:size(BGSeg,2)
    BGP1(:,i) = PG.get('1pos',BGSeg(1,i));
    BGP2(:,i) = PG.get('1pos',BGSeg(2,i));
    BGPrel(:,i) = BGP2(:,i)-BGP1(:,i);
    BGPhi(i) = atan2(BGPrel(2,i),BGPrel(1,i));
end
BGSigma = (diag(BGPrel.'*BGPrel).^0.5).';
[BGDepth,BGSigma,EscapeFeatures] = DSEscape(PG,BGSeg,BGSigma,BGPhi);
timeDS = toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure;
hold on;
plot(BGS,BGDepth);
xlabel('basket grasp curve parameter $u$','interpreter','latex');
ylabel('$double-support \Delta U$','interpreter','latex');
% plot(BGS,EscapeFeatures);


%% calculating the SS depth graph and updating the bg depth graph
tic
BGDepth2 = BGDepth;
SSDepth = ones(size(BGDepth))*inf;
SSEscape_Features = zeros(size(EscapeFeatures));
[bg,sigmabg,~] = bgheight(PG,BGSegEnds);
for i=1:2*PG.nv
    featurepair = (i<=PG.nv)*[0,i]+(i>PG.nv)*[i-PG.nv,0];
    [ppv,sigma1] = depthperfeature2(PG,BGSigma,BGPhi,featurepair);
    interpbg = linspace(bg(2,1),bg(2,2),length(sigma1));
    tempSSDepth = ppv-interpbg;
    [SSDepth,Ind] = min([SSDepth;max([zeros(size(ppv));tempSSDepth],[],1)],[],1);
    SSEscape_Features(:,Ind) = ones(1,nnz(Ind)).*featurepair.';
    if all(tempSSDepth>0) && any(tempSSDepth<BGDepth2)
        [BGDepth2,Index] = min([BGDepth2;SSDepth],[],1);
    end

    %         if BGDepth2(1)<1
    %             disp('here')
    %         end
end
% need to make a check before SS computation of subsections of BGS where
% BGDepth = inf, to remove them from the computation (caged grasps have
% infinite depth)
for i=1:length(BGDepth) %depth for caged grasps is still infinite (replacing a real check for now)
    if BGDepth(i) == inf
        BGDepth2(i) = inf;
        SSDepth(i) = inf;
    end
end
timeSS = toc
%% Comparing with grasp search computation of final escape depth (Slow!)
if 0
    N = 10;
    BGs1 = linspace(BGSegEnds(1,1),BGSegEnds(1,2),N);
    BGs2 = linspace(BGSegEnds(2,1),BGSegEnds(2,2),N);
    BGslist = [BGs1;BGs2];
    BGdepthlist = zeros(1,size(BGslist,2));
    sig = zeros(1,size(BGslist,2));
    pplist = zeros(size(BGslist,2),7);
    phi = zeros(1,size(BGslist,2));
    tic
    for i=1:N
        BGSv = BGslist(:,i);
        [BGdepthlist(i),sig(i),phi(i),pp] = BasketDepth(PG,BGSv);
        pplist(i,1:length(pp)) = pp;
    end
    BGSn = linspace(0,1,N);
end
%%
figure;
hold on;

plot(BGS,BGDepth2,'g','lineWidth',3,'DisplayName','Local escape');
% plot(BGSn,BGdepthlist,'+k','MarkerSize',8,'lineWidth',3,'DisplayName','Graph search');
plot(BGS,BGDepth,'r--','lineWidth',2,'DisplayName','Double-support');
plot(BGS,SSDepth,'b--','lineWidth',2,'DisplayName','Single-support');
% title('Basket Depth for Basket grasp segment - Cup','interpreter','latex','fontSize',14);
xlabel('u','interpreter','latex');
ylabel('$\Delta U$','interpreter','latex');
legend()
hold off

%% Saving data for basket grasp curve segment, in case of multi-segment curve
% save('SegGraphData.mat','BGSigma','BGDepth2','BGDepth','SSDepth','BGSegEnds')
%% Plotting the graph of multi-segment basket grasp curve
if 0
figure
plotmultisegment('NewSeg1GraphData.mat','NewSeg4GraphData.mat','NewSeg3GraphData.mat')
end
%% Some extra code snippets for analyzing interesting points along the basket grasp curve
hold on
interestingList = [0.1952,0.6817];
for i=1:length(interestingList)
plot(interestingList(i)*ones(2,1),[0,0.55],'k--')
end

%% Find the analytical segments of the basket grasps

% tic
% Contours = equisigmacontourSym(PG,BGSegEnds);
% toc
% disp('success!!')


%% sample the grasp curve according to a normal distribution to find mean
samp_bgu = 0.5*randn(5000000,1);
sum(abs(samp_bgu)<min(abs(min(BGS)),abs(max(BGS))))
samp_bgu(abs(samp_bgu)>=1)=[];
samp_bgu = samp_bgu/2+0.5;
samp_depth = interp1(BGS,BGDepth,samp_bgu);
figure
plot(samp_bgu,samp_depth,'.')
figure
histogram(samp_bgu);
max(samp_depth)/mean(samp_depth)
%%
target = find(BGS>=0.801,1,'first');
targetmin = find(BGS>=0.801-0.18,1,'first');
targetmax = find(BGS>=0.801+0.02,1,'first');
depth_imp = min(BGDepth(targetmin:targetmax))
% depth_imp = max(BGDepth);
relative_security1 = depth_imp/BGDepth(end/2)
relative_security2 = depth_imp/mean(samp_depth)
relative_security3 = depth_imp/mean(BGDepth)




%% draw the gravitational contact space for a point from the contours
% n = 10;
% Contour = Contours{n}{1};
% segment = Contour{1};
% sinit = [17.07,25.58];%segment{end-1};
% n =round((interestingList(2)+0.1)*reso,0);
% n =round(0.5*reso,0);
% 
% sinit = BGSeg(:,n);
sinit = BGslist(:,30);
if 1
[f1,f2] = PG.get('2pos',sinit(1),sinit(2));
basepos = [f1,f2];
sig = norm(f2-f1);
% smaxv = [smaxv1(:,n),smaxv2(:,n)];
% thetamaxv = [thetamaxv1(n),thetamaxv2(n)];
% hmaxv = [hmaxv1(n),hmaxv2(n)];
% hmaxind = find(hmax == min(hmax));
% smax = smaxv(:,hmaxind);
% thetamax = thetamaxv(hmaxind);
% hmax = hmaxv(hmaxind);
[cont_original] = PG.GetSigmaContours(Sigma,sig);


cont = PG.CleanContour(cont_original,basepos);

[ds_max,ds_min,ds_virtual] = PG.DSNodes(cont);
[ss_max,ss_min,ss_saddle] = PG.SSNodes();
nodes{1} = ds_max;
nodes{2} = ds_min;
nodes{3} = ds_virtual;
nodes{4} = ss_max;
nodes{5} = ss_min;
nodes{6} = ss_saddle;

nodes = check_SS_for_DS(PG,nodes,f1,f2); %function checks SS_nodes for penetration (of the other finger),
% each node appears once for each finger it is relevant for, with the relevant finger index at the end 
% save([pwd '/' subFolderName '/PreGraph.mat'],'PG','nodes','cont','starting_node','f1','f2');
%% draw the contact space
loadfromfile = 0;
if loadfromfile == 1
    for finger = 1:2
        if exist([pwd '/' subFolderName '/Nodes s' num2str(finger) '.fig'],'file')
            graphfig{finger} = openfig([pwd '/' subFolderName '/Nodes s' num2str(finger) '.fig']); 
%             graphfig{finger} = openfig([pwd '/' subFolderName '/Polygon Escape Graph' num2str(finger) '.fig']); 
        end
    end
end
if ~exist('graphfig','var')
    graphfig = cell(1,2);
end

for finger = 1:2
    if isempty(graphfig{finger})||~isgraphics(graphfig{finger})
            otherFingerRelative = basepos(:,3-finger)-basepos(:,finger);
            contactHeight = basepos(2,finger);
            graphfig{finger} = plotGraphNodes(PG,cont,finger,nodes,otherFingerRelative,contactHeight);
            hold on
%             plot(smax(finger),thetamax,'k+','markerSize',14,'lineWidth',6)
            hold off
    end
    set(graphfig{finger},'Position',[960*2+(finger-1)*960,41,960,920])
    set(gcf,'renderer','Painters')
end

end
%% plot the BG and escape grasps together for specific sigma values
% sigmaval = 6.45;
% BGSval = interp1(sig,BGslist.',sigmaval);
if 1
BGSval = BGslist(:,96);
[BGdepthval,sigval,phival,ppval] = BasketLocalDepth(PG,BGSval);
[BGdepthval,sigval] = BasketLocalDepthPath(PG,BGSval,[]);

PG.drawPolygonClear()
[f1pos,f2pos] = PG.get('2Pos',BGSval(1),BGSval(2));
hold on
if length(ppval)>4
drawPolygonNodePrepare(PG,ppval(1),ppval(2),ppval(4),f1pos,f2pos)
else
    s1 = [];s2 = [];
    if ppval(end)==1
        s1 = ppval(1);
    else
        s2 = ppval(1);
    end
    drawPolygonNodePrepare(PG,s1,s2,ppval(3),f1pos,f2pos)
end

hold on
end

%%
% figure;hold on;
% plot(1:N,BGdepthlist);
% plot(1:N,sig);hold off;
figure;
hold on
p1 = plot(sig(1:N),BGdepthlist(1:N),'.k','markersize',10,'DisplayName','Algorithm');
% plot(sig,phi,'.r','markersize',10);
p2 = plot(sig(switchlist),BGdepthlist(switchlist),'+r','markersize',10,'DisplayName','Escape-feature change');
xlabel('$\sigma$','interpreter','latex','fontsize',25)
ylabel('$\Delta U$','interpreter','latex','fontsize',25)
title('Basket depth vs. inter-finger distance','interpreter','latex','fontsize',20)

labelflag = 0;
% for i=1:numel(depthcurves)
%     depthc = depthcurves{i};
%     if ~isempty(depthc)
%         if labelflag == 0
%             p4 = plot(depthc(2,:),depthc(1,:),'b','DisplayName','Analytical depth');
%             labelflag = 1;
%         else
%             plot(depthc(2,:),depthc(1,:),'b')
%         end
%     end
% end
labelflag = 0;
for i=1:numel(depthcurves2)
    depthc = depthcurves2{i};
    if ~isempty(depthc)
        if labelflag == 0
            p3 = plot(depthc(2,:),depthc(1,:),'g','DisplayName','other feature pairs');
            labelflag = 1;
        else
        plot(depthc(2,:),depthc(1,:),'g')
        end
    end
end
% legend([p1,p2,p3,p4])
legend([p1,p2,p3])
% legend('Algorithm','Escape-feature change','other feature pairs','Analytical depth')
hold off



%% check specific sigma value
% i = 2;
% BGS = BGslist(:,i);
% [BGdepthlist(i),sig(i)] = BasketDepthPath(PG,BGS,4);
%% Draw caging regions
% cagedpoly = PG.drawPolygon();
% hold on;
% for i=2:1:N
%     BGS = BGslist(:,i);
%     dep = BGdepthlist(i);
%     sigi = sig(i);
%     col = 'g';
%     drawCaging(PG,BGS,dep,col);
% end
% for i=46:1:50
%     BGS = BGslist(:,i);
%     dep = BGdepthlist(i);
%     sigi = sig(i);
%     col = 'r';
%     drawCaging(PG,BGS,dep,col);
% end
%%
% cagedpoly2 = PG.drawPolygon();
% hold on;
% for i=10:1:N
%     BGS = BGslist(:,i);
%     dep = BGdepthlist(i)-0.2;
%     sigi = sig(i);
%     col = 'g';
%     drawCaging(PG,BGS,dep,col);
% end
function [bg,sigma,phi] = bgheight(PG,BGslist)
% assuming a single edge for each finger, find their parameters:
edgenum1 = PG.get('edgeNum',BGslist(1,1));
edgenum2 = PG.get('edgeNum',BGslist(2,1));
t1 = PG.tangent(:,edgenum1);
t2 = PG.tangent(:,edgenum2);
t1 = t1/norm(t1);
t2 = t2/norm(t2);
start1 = PG.get('1Pos',BGslist(1,1));
start2 = PG.get('1Pos',BGslist(2,1));

deltaS = BGslist-BGslist(:,1); %get each point as a delta-s relative to the start
bg = -(start1+deltaS(1,:).*t1); %find the COM relative to the first finger
p1 = start1+deltaS(1,:).*t1;
p2 = start2+deltaS(2,:).*t2;
distv = (p2-p1);
sigma = sqrt(distv(1,:).^2+distv(2,:).^2);
phi = atan2(distv(2,:),distv(1,:));
end


function vertex = findvertex(PG,s)
vertex1 = find(PG.S(PG.VL)==s(1));
vertex2 = find(PG.S(PG.VL)==s(2));
if ~isempty(vertex1)
    vertex = [1,vertex1];
    return
end
if ~isempty(vertex2)
    vertex = [2,vertex2];
    return
end
vertex1 = find(abs(PG.S(PG.VL)-s(1))==min(abs((PG.S(PG.VL)-s(1)))));
vertex2 = find(abs(PG.S(PG.VL)-s(2))==min(abs((PG.S(PG.VL)-s(2)))));
if min(abs((PG.S(PG.VL)-s(2))))<min(abs((PG.S(PG.VL)-s(1))))
    vertex = [2,vertex2];
else
    vertex = [1,vertex1];
end
if min(min(abs((PG.S(PG.VL)-s(2)))),min(abs((PG.S(PG.VL)-s(1)))))>1E-6
    disp('not on rectangle edge');
end
end

