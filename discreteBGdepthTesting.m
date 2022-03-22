clearvars
close all
%% Define the polygon
% %% Define the polygon
clc
object = 4;
symmetry = 0;
switch object
    case 0 % pentagon
        theta = 0;%-pi()/2;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
        PG = Polygon(R*[3,-4;3,4;-1.5,5;-3,0;-1.5,-5].'+[0;3],[0;3]);
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
        object = [-2,-1;2,-1;3,3;1,0;0,1;-1,0;-4,1].';
        com = [0,0].';
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
%          object = round(object./20,1);
%          com = round(com./20,0);
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
end

polyDrawing = PG.drawPolygon();
hold off
%% Find equilibrium grasps

[PG,S,X,VL] = PG.findBdyVariable(res);
[s1,s2,~] = PG.Eqcheck();

%% Check grasps for minimum/maximum
[BG,NBG] = PG.StatusSeparate(s1,s2,X);
%% Full contact space figure
contspace = figure;
hold on
c_layers=20;
Sigma=inter_finger_distance(X,X);
contour(S,S,Sigma,c_layers)
xlabel('$s_1$','interpreter','latex','fontsize',18)
ylabel('$s_2$','interpreter','latex','fontsize',18)
colorbar
col = [0 1 0; 0 0 0; 1 0 0];
plot(BG(1,:),BG(2,:),'.','MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(1,:),'markersize',8)
plot(NBG(1,:),NBG(2,:),'.','MarkerEdgeColor',col(3,:),'MarkerFaceColor',col(3,:),'markersize',8)
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
plot(BGSegEnds(1,:),BGSegEnds(2,:),'m+','markerSize',14,'lineWidth',5)
hold off
%%
drawPolygon2fGraspRange(PG,BGSegEnds)

%% Find basket grasp depth at discrete points along the segments:
N = 20;
BGs1 = linspace(BGSegEnds(1,1),BGSegEnds(1,2),N);
BGs2 = linspace(BGSegEnds(2,1),BGSegEnds(2,2),N);
BGslist = [BGs1;BGs2];
% BGslist(:,[1,N+2]) = [];01
BGdepthlist = zeros(1,size(BGslist,2));
sig = zeros(1,size(BGslist,2));
pplist = zeros(size(BGslist,2),7);
phi = zeros(1,size(BGslist,2));
tic
for i=1:N
    BGS = BGslist(:,i);
    [BGdepthlist(i),sig(i),phi(i),pp] = BasketDepth(PG,BGS);
    pplist(i,1:length(pp)) = pp;
end
toc

%% plot the BG and escape grasps together for specific sigma values
sigmaval = 6.45;
BGSval = interp1(sig,BGslist.',sigmaval);
[BGdepthval,sigval,phival,ppval] = BasketDepth(PG,BGSval.');
% [BGdepthval,sigval] = BasketDepthPath(PG,BGSval.',[]);

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
%% find the analytical depth of the feature-segment
% recieve the puncture point from the basket depth algorithm
% identify switches between different feature types (V-E,E-E,E)
switchlist = 1; 
featurelist = zeros(size(BGslist,2),2); % for each i, both fingers contact feature is defined, negative for vertices, positive for edges, 0 for non-contact
% example - (0,1) means a SS with O2 supporting the first edge. 
% (-2,1) means DS, with O1 supporting the second V and O2 the first E.
for i=1:N
    pp = pplist(i,:);
    numpar = find(pp~=0,1,'last');
    pp = pp(1:numpar);
    %     BGS = BGslist(:,i);
    feature = [0,0];
    changed = 0;
    switch length(pp)

        case 4 % E / V
            vertnum = find(PG.S(PG.VL)==pp(1));
            if ~isempty(vertnum) % V
                feature(pp(4)) = -vertnum;
                if length(symmetry)>1
                    fchange = symmetry(3,-feature(pp(4))); % change only if necessary
                    if fchange~=feature(pp(4))
                        feature(3-pp(4)) = fchange;
                        feature(pp(4)) = 0;
                    end
                end
            else % E
                edgenum = PG.get('edgeNum',pp(1));
                feature(pp(4)) = edgenum;
                if length(symmetry)>1 % correct the features according to symmetry
                    fchange = symmetry(1,feature(pp(4))); % change only if necessary
                    if fchange~=feature(pp(4))
                        feature(3-pp(4)) = fchange;
                        feature(pp(4)) = 0;
                    end
                end
            end
        case {6,7} % E-V / V-E / E-E / V-V
            % first finger:
            feature(1) = PG.get('edgeNum',pp(1));
            if pp(1) == PG.S(PG.VL(feature(1)))
                feature(1) = -feature(1);
                if length(symmetry)>1
                    fchange = symmetry(3,-feature(1)); % change only if necessary
                    if fchange~=feature(1)
                        changed = 1;
                        feature(2) = fchange;
                    end
                end
            else
                if length(symmetry)>1
                    fchange = symmetry(1,feature(1)); % change only if necessary
                    if fchange~=feature(1)
                        changed = 1;
                        feature(2) = fchange;
                    end
                end
            end
            
            % second finger:
            f2 = PG.get('edgeNum',pp(2));
            if pp(2) == PG.S(PG.VL(f2))
                f2 = -f2;
                if changed
                    feature(1) = symmetry(4,-f2); % change to symmetrical vertex
                end
            else
                if changed
                    feature(1) = symmetry(2,f2); % change to symmetrical edge
                end
            end
            if feature(2) == 0 % in case no symm. changes were made
                feature(2) = f2;
            end
            
    end
    featurelist(i,:) = feature;
    if i>1
        if any(featurelist(i,:)~=featurelist(i-1,:))
            switchlist = [switchlist,i]; %#ok<AGROW>
        end   
%         end
    end
end


% construct the curves for each feature type
tic
[bg,sigmabg,phibg] = bgheight(PG,BGslist);
if sigmabg(end)<sigmabg(1)
    sigmabg = sigmabg(end:-1:1);
    bg = bg(:,end:-1:1);
    phibg = phibg(end:-1:1);
end
depthcurves = cell(size(switchlist));
for i=1:length(switchlist)
    featurepair = featurelist(switchlist(i),:);
    reso=100;
    [ppv,sigma] = depthperfeature(PG,sigmabg,phibg,featurepair,reso);
    interpbg = interpn(sigmabg,bg(2,:),sigma);
    depthcurves{i} = [ppv-interpbg;sigma];
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% %% find each depth curve for all possible escape feature
% tic
% possiblefeaturelist = [];
% for i=1:PG.nv
%     possiblefeaturelist = [possiblefeaturelist;i 0;0 i;-i 0;0 -i];
%     for j=1:PG.nv
%         addition = [];
%         if i~=j
%             addition = [-i,j;j,-i;i,j];
%         end
%         possiblefeaturelist = [possiblefeaturelist;addition]; 
%     end
% end
% depthcurves2 = cell(size(possiblefeaturelist,2),1);
% for i=1:size(possiblefeaturelist,1)
%     featurepair = possiblefeaturelist(i,:);
%     reso=100;
%     [ppv,sigma] = depthperfeature(PG,sigmabg,phibg,featurepair,reso);
%     interpbg = interpn(sigmabg,bg(2,:),sigma);
%     depthcurves2{i} = [ppv-interpbg;sigma];
% end
% toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find each depth curve iteratively, starting from the basket grasp feature pair:
tic
possiblefeaturelist = [];
% positive direction:

possiblefeaturelist = [-7 2;7 -2;7 1;-7 1];

depthcurves2 = cell(size(possiblefeaturelist,2),1);
for i=1:size(possiblefeaturelist,1)
    featurepair = possiblefeaturelist(i,:);
    reso=100;
    [ppv,sigma] = depthperfeature(PG,sigmabg,phibg,featurepair,reso);
    interpbg = interpn(sigmabg,bg(2,:),sigma);
    depthcurves2{i} = [ppv-interpbg;sigma];
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function [pp,sigma] = depthperfeature(PG,sigmaseg,phiseg,featurepair,N)
nonz = find(featurepair~=0);
switch length(nonz)
    case 1
        if featurepair(nonz)>0
            [pp,sigma] = depthperfeatureE(PG,sigmaseg,phiseg,featurepair,N);
        else
            [pp,sigma] = depthperfeatureV(PG,sigmaseg,phiseg,featurepair,N);
        end
    case 2
        if featurepair(1)*featurepair(2)<0 %one finger supports a vertex
            [pp,sigma] = depthperfeatureEV(PG,sigmaseg,phiseg,featurepair,N);
        else
            [pp,sigma] = depthperfeatureEE(PG,sigmaseg,phiseg,featurepair,N);
        end
end
end

function [pp,sigma] = depthperfeatureEV(PG,sigmaseg,phiseg,featurepair,N)
% function recieves a segment of sigma and phi values, a feature pair
% (negative for vertex, positive for edge), a resolution N for the curve
% function returns a vector of puncture point heights, relative to the
% first finger, over a possible range of sigma
vertfinger = find(featurepair<0);
vertnum = -featurepair(vertfinger);
edgefinger = 3-vertfinger;
edgenum = featurepair(edgefinger);
% find sigma limits for feature pair:
vertic = [PG.vertex PG.vertex(:,1)];
dist1 = norm(vertic(:,vertnum)-vertic(:,edgenum)); %distance between vertices
dist2 = norm(vertic(:,vertnum)-vertic(:,edgenum+1));
a = vertic(:,edgenum+1)-vertic(:,edgenum);
b = vertic(:,edgenum)-vertic(:,vertnum);
dist3 = abs(a(1)*b(2)-a(2)*b(1))/norm(a); % distance between vertex and edge
minsig = min([dist1,dist2,dist3]);
maxsig = max([dist1,dist2,dist3]);
maxsig = min(maxsig,max(sigmaseg));
minsig = max(minsig,min(sigmaseg));
N = round(N*(maxsig-minsig));
sigma = linspace(minsig,maxsig,N); % segment of relevant sigmas
phi = interpn(sigmaseg,phiseg,sigma);
if vertfinger == 2 %switch the angles to work around finger 2, fix in the end
    phi = phi+pi();
end
% find puncture point com position
vj = PG.vertex(:,vertnum);
vk = PG.vertex(:,edgenum);
tk = PG.tangent(:,edgenum);
% nk = PG.normal(:,edgenum);
% calculate the range of angles in which the normal from the supported vertex must be
nj2 = PG.normal(:,vertnum);
if vertnum>1
nj1 = PG.normal(:,vertnum-1);
else
    nj1 = PG.normal(:,end);
end
vertangle1 = atan2(nj1(2),nj1(1));
vertangle2 = atan2(nj2(2),nj2(1));
% check if vertex is convex
J = [0 1;-1 0];
tj2 = J*nj2;
if nj1.'*tj2>0 %convex
    posbaseangle = vertangle1;
    posrangeangle = wrapToPi(vertangle2-vertangle1);
else
    posbaseangle = vertangle2;
    posrangeangle = wrapToPi(vertangle1-vertangle2);
end
pp = [];


% R = [cos(theta) -sin(theta);sin(theta) cos(theta)];

p2 = [sigma.*cos(phi);sigma.*sin(phi)];

dp1 = p2(1,:);
dp2 = p2(2,:);
vj1 = vj(1);
vj2 = vj(2);
vk1 = vk(1);
vk2 = vk(2);
tk1 = tk(1);
tk2 = tk(2);
%calculate (in parallel) the escape stance
cyv = [(dp2.*vj1.^2 + dp2.*vj2.^2 + dp1.*vj1.*vk2 - dp1.*vj2.*vk1 - dp2.*vj1.*vk1 - dp2.*vj2.*vk2)./(dp1.^2 + dp2.^2) + ((dp1.*tk1.*vj2 - dp1.*tk2.*vj1 + dp2.*tk1.*vj1 + dp2.*tk2.*vj2).*(tk1.*vk1 - tk2.*vj2 - tk1.*vj1 + tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2)))./((tk1.^2 + tk2.^2).*(dp1.^2 + dp2.^2));...
    (dp2.*vj1.^2 + dp2.*vj2.^2 + dp1.*vj1.*vk2 - dp1.*vj2.*vk1 - dp2.*vj1.*vk1 - dp2.*vj2.*vk2)./(dp1.^2 + dp2.^2) - ((dp1.*tk1.*vj2 - dp1.*tk2.*vj1 + dp2.*tk1.*vj1 + dp2.*tk2.*vj2).*(tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2)))./((tk1.^2 + tk2.^2).*(dp1.^2 + dp2.^2))];
%     thetav = [ 2.*atan((dp1.*tk1 + dp2.*tk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1));...
%             2.*atan((dp1.*tk1 + dp2.*tk2 - (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))];
av = [-(tk1.*vk1 - tk2.*vj2 - tk1.*vj1 + tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(tk1.^2 + tk2.^2);...
    (tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(tk1.^2 + tk2.^2)];
thetav = [2.*atan((dp1.*tk1 + dp2.*tk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1));...
 2.*atan((dp1.*tk1 + dp2.*tk2 - (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))];

% for eq. condition, calculate the normal intersection point, and its angle
% to support 1
% betav1 = (tk1.*vj1 + tk2.*vj2 - (dp1^2.*tk1^2 + dp1^2.*tk2^2 + dp2^2.*tk1^2 + dp2^2.*tk2^2 - tk1^2.*vj2^2 + 2.*tk1^2.*vj2.*vk2 - tk1^2.*vk2^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2^2.*vj1^2 + 2.*tk2^2.*vj1.*vk1 - tk2^2.*vk1^2)^(1./2))./(tk2.*cos(2.*atan((dp1.*tk1 + dp2.*tk2 + (dp1^2.*tk1^2 + dp1^2.*tk2^2 + dp2^2.*tk1^2 + dp2^2.*tk2^2 - tk1^2.*vj2^2 + 2.*tk1^2.*vj2.*vk2 - tk1^2.*vk2^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2^2.*vj1^2 + 2.*tk2^2.*vj1.*vk1 - tk2^2.*vk1^2)^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))) + tk1.*sin(2.*atan((dp1.*tk1 + dp2.*tk2 + (dp1^2.*tk1^2 + dp1^2.*tk2^2 + dp2^2.*tk1^2 + dp2^2.*tk2^2 - tk1^2.*vj2^2 + 2.*tk1^2.*vj2.*vk2 - tk1^2.*vk2^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2^2.*vj1^2 + 2.*tk2^2.*vj1.*vk1 - tk2^2.*vk1^2)^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))));
angleFromSupp1 = angle(vj2.*sin(2.*atan((tk1.*vk1 - tk2.*vj2 - tk1.*vj1 + tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1) + (dp1.*tk1 + dp2.*tk2 + tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2)./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))) - vj2.*cos(2.*atan((tk1.*vk1 - tk2.*vj2 - tk1.*vj1 + tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1) + (dp1.*tk1 + dp2.*tk2 + tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2)./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))).*1i - vj1.*sin(2.*atan((tk1.*vk1 - tk2.*vj2 - tk1.*vj1 + tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1) + (dp1.*tk1 + dp2.*tk2 + tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2)./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))).*1i - vj1.*cos(2.*atan((tk1.*vk1 - tk2.*vj2 - tk1.*vj1 + tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1) + (dp1.*tk1 + dp2.*tk2 + tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2)./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))) + ((tk1.*vj1 + tk2.*vj2 - (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2)).*1i)./(tk2.*cos(2.*atan((dp1.*tk1 + dp2.*tk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))) + tk1.*sin(2.*atan((dp1.*tk1 + dp2.*tk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1)))));
angleFromSupp1(~isreal(thetav(1,:))) = NaN;
try
relangleSupp1 = wrapToPi(angleFromSupp1-(posbaseangle+thetav(1,:))); %%%%%%%%NEED to fix this, add the theta angle there (with apropriate sign)
catch
    disp('imaginary error');
end
condangle11 = relangleSupp1-posrangeangle<=0 & relangleSupp1>=0;
relanglesupp1neg = wrapToPi(angleFromSupp1-(posbaseangle+thetav(1,:)+pi()));
condangle12 = relanglesupp1neg-posrangeangle<=0 & relanglesupp1neg>=0;
condangle1 = condangle11 | condangle12;
% % njeff1 = [cos(angleFromSupp1);sin(angleFromSupp1(1,:))];
% % lamda1 = [njeff1 nk]\[0;1];%%%%%%%%NEED to fix this, add the theta angle there (with apropriate sign)
lamda1 =  [(tk2.*cos(thetav(1,:)) + tk1.*sin(thetav(1,:)))./(tk1.*cos(angleFromSupp1).*cos(thetav(1,:)) - tk2.*cos(angleFromSupp1).*sin(thetav(1,:)) + tk2.*sin(angleFromSupp1).*cos(thetav(1,:)) + tk1.*sin(angleFromSupp1).*sin(thetav(1,:)));...
                cos(angleFromSupp1)./(tk1.*cos(angleFromSupp1).*cos(thetav(1,:)) - tk2.*cos(angleFromSupp1).*sin(thetav(1,:)) + tk2.*sin(angleFromSupp1).*cos(thetav(1,:)) + tk1.*sin(angleFromSupp1).*sin(thetav(1,:)))];
condlamda1 = all(lamda1>0);
% % njeff2 = [cos(angleFromSupp1(2,:));sin(angleFromSupp1(2,:))];
% % lamda2 = [njeff2 nk]\[0;1];
% 
% betav2 = (tk1.*vj1 + tk2.*vj2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(tk2.*cos(2.*atan((dp1.*tk1 + dp2.*tk2 - (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))) + tk1.*sin(2.*atan((dp1.*tk1 + dp2.*tk2 - (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))));
angleFromSupp2 = angle(- vj1.*cos(2.*atan((tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1) - (dp1.*tk1 + dp2.*tk2 + tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2)./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))) - vj2.*cos(2.*atan((tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1) - (dp1.*tk1 + dp2.*tk2 + tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2)./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))).*1i + vj1.*sin(2.*atan((tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1) - (dp1.*tk1 + dp2.*tk2 + tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2)./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))).*1i - vj2.*sin(2.*atan((tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1) - (dp1.*tk1 + dp2.*tk2 + tk1.*vj1 + tk2.*vj2 - tk1.*vk1 - tk2.*vk2)./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))) + ((tk1.*vj1 + tk2.*vj2 + (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2)).*1i)./(tk2.*cos(2.*atan((dp1.*tk1 + dp2.*tk2 - (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1))) + tk1.*sin(2.*atan((dp1.*tk1 + dp2.*tk2 - (dp1.^2.*tk1.^2 + dp1.^2.*tk2.^2 + dp2.^2.*tk1.^2 + dp2.^2.*tk2.^2 - tk1.^2.*vj2.^2 + 2.*tk1.^2.*vj2.*vk2 - tk1.^2.*vk2.^2 + 2.*tk1.*tk2.*vj1.*vj2 - 2.*tk1.*tk2.*vj1.*vk2 - 2.*tk1.*tk2.*vj2.*vk1 + 2.*tk1.*tk2.*vk1.*vk2 - tk2.^2.*vj1.^2 + 2.*tk2.^2.*vj1.*vk1 - tk2.^2.*vk1.^2).^(1./2))./(dp1.*tk2 - dp2.*tk1 + tk1.*vj2 - tk2.*vj1 - tk1.*vk2 + tk2.*vk1)))));
angleFromSupp2(~isreal(thetav(2,:))) = NaN;
relangleSupp2 = wrapToPi(angleFromSupp2-(posbaseangle+thetav(2,:)));
condangle21 = relangleSupp2-posrangeangle<=0 & relangleSupp2>=0;
relanglesupp2neg = wrapToPi(angleFromSupp2-(posbaseangle+thetav(2,:)+pi()));
condangle22 = relanglesupp2neg-posrangeangle<=0 & relanglesupp2neg>=0;
condangle2 = condangle21 | condangle22;
% % njeff2 = [cos(angleFromSupp2);sin(angleFromSupp2)];
% % lamda2 = [njeff2 nk]\[0;1];
lamda2 =  [(tk2.*cos(thetav(2,:)) + tk1.*sin(thetav(2,:)))./(tk1.*cos(angleFromSupp2).*cos(thetav(2,:)) - tk2.*cos(angleFromSupp2).*sin(thetav(2,:)) + tk2.*sin(angleFromSupp2).*cos(thetav(2,:)) + tk1.*sin(angleFromSupp2).*sin(thetav(2,:)));...
                cos(angleFromSupp2)./(tk1.*cos(angleFromSupp2).*cos(thetav(2,:)) - tk2.*cos(angleFromSupp2).*sin(thetav(2,:)) + tk2.*sin(angleFromSupp2).*cos(thetav(2,:)) + tk1.*sin(angleFromSupp2).*sin(thetav(2,:)))];
condlamda2 = all(lamda2>0);
cond3 = [condlamda1;condlamda2];

cond2 = [condangle1;condangle2]; %%% added Equilibrium conditions, ADD afterwards saddle-condition
cond1 = (av<=1&av>=0&imag(av)==0);
cond = cond1 & cond2 & cond3;
% can derive the determinant of the grasp matrix, and calculate it
% vectorically.
solutionSt = find(any(cond),1,'first');
solutionEn = find(all(~cond(:,solutionSt:end)),1,'first');
if isempty(solutionEn)
    solutionEn = size(cond,2);
else
    solutionEn = solutionEn+solutionSt-2;
end
for i=1:(solutionEn-solutionSt+1)
    try
    pp(i) = cyv(cond(:,solutionSt+i-1),solutionSt+i-1);
    catch me
        if diff(av(:,solutionSt+i-1))==0
            pp(i) = cyv(1,solutionSt+i-1);
        else
                    disp(me)
        end

    end
end
%     pp = cyv(solution);
if vertfinger == 2 %(need to correct for finger movement)
    pp = pp-p2(2,solutionSt:solutionEn); %add the second finger's height
end
sigma = sigma(solutionSt:solutionEn);
end

function [pp,sigma] = depthperfeatureEE(PG,sigmaseg,phiseg,featurepair,N)
% function recieves a segment of sigma and phi values, a feature pair
% (of two edge numbers), a resolution N for the curve
% function returns a vector of puncture point heights, relative to the
% first finger, over a possible range of sigma
pp = [];
sigma = [];
edge1 = featurepair(1);
edge2 = featurepair(1);
dels1 = PG.S(PG.VL(edge1+1))-PG.S(PG.VL(edge1));
dels2 = PG.S(PG.VL(edge2+1))-PG.S(PG.VL(edge2));
        nj = PG.normal(:,edge1);
        nj = nj/norm(nj);
        nk = PG.normal(:,edge2);
        nk = nk/norm(nk);
        epsi = 1E-7;
        J = [0 -1;1 0];
        if abs(nj.'*J*nk)<epsi
            return
        end
        vj = PG.vertex(:,edge1);
        vk = PG.vertex(:,edge2);
        %
        Sv = linspace(Sstart,Send,N);
        sol1 = zeros(3,size(Sv,2));
        sol2 = zeros(3,size(Sv,2));
        sol3 = zeros(3,size(Sv,2));
        sol4 = zeros(3,size(Sv,2));
        dp = zeros(2,size(Sv,2));
        for i=1:N
            dp(:,i) = sigmaseg.*[cos(phiseg);sin(phiseg)];
            [sol1(:,i),sol2(:,i),sol3(:,i),sol4(:,i)] = EEheight(dp(:,i),nj,nk,vj,vk,dels1,dels2);
        end

        
        if any(~isnan(sol1(3,:)))
            curveSt = find(~isnan(sol1(3,:)),1,'first');
            curveEn = curveSt + find(isnan(sol1(3,curveSt+1:end)),1,'first')-1;
            curve = [sol1(:,curveSt:curveEn);sigmaseg(curveSt:curveEn)];
            if ~isempty(find(~isnan(sol1(3,curveEn+1:end)),1))
                disp('Discontinuity!')
            end
        end
        if any(~isnan(sol2(3,:)))
            curveSt = find(~isnan(sol2(3,:)),1,'first');
            curveEn = curveSt + find(isnan(sol2(3,curveSt+1:end)),1,'first')-1;
            curve = [sol2(:,curveSt:curveEn);sigmaseg(curveSt:curveEn)];
            if ~isempty(find(~isnan(sol2(3,curveEn+1:end)),1))
                disp('Discontinuity!')
            end
        end
        if any(~isnan(sol3(3,:)))
            curveSt = find(~isnan(sol3(3,:)),1,'first');
            curveEn = curveSt + find(isnan(sol3(3,curveSt+1:end)),1,'first')-1;
            curve = [sol3(:,curveSt:curveEn);sigmaseg(curveSt:curveEn)];
            if ~isempty(find(~isnan(sol3(3,curveEn+1:end)),1))
                disp('Discontinuity!')
            end
        end
        if any(~isnan(sol4(3,:)))
            curveSt = find(~isnan(sol4(3,:)),1,'first');
            curveEn = curveSt + find(isnan(sol4(3,curveSt+1:end)),1,'first')-1;
            curve = [sol4(:,curveSt:curveEn);sigmaseg(curveSt:curveEn)];
            if ~isempty(find(~isnan(sol4(3,curveEn+1:end)),1))
                disp('Discontinuity!')
            end
        end
        pp = curve(2,:);
        sigma = curve(end,:);
end

function [pp,sigma] = depthperfeatureV(PG,sigmaseg,phiseg,featurepair,N)

vertfinger = find(featurepair<0);
vertnum = -featurepair(vertfinger);
dist = norm(PG.vertex(:,vertnum));
sigma = linspace(sigmaseg(1),sigmaseg(end),abs(round(N*(sigmaseg(1)-sigmaseg(end)))));
phi = interpn(sigmaseg,phiseg,sigma);
pp = ones(size(sigma))*dist;
if vertfinger == 2
    pp = pp + sigma.*sin(phi);
end
end

function [pp,sigma] = depthperfeatureE(PG,sigmaseg,phiseg,featurepair,N)

edgefinger = find(featurepair>0);
edgenum = featurepair(edgefinger);
vertic = [PG.vertex PG.vertex(:,1)];
a = vertic(:,edgenum+1)-vertic(:,edgenum);
b = vertic(:,edgenum);
dist = abs(a(1)*b(2)-a(2)*b(1))/norm(a); % distance between vertex and edge
sigma = linspace(sigmaseg(1),sigmaseg(end),abs(round(N*(sigmaseg(1)-sigmaseg(end)))));
phi = interpn(sigmaseg,phiseg,sigma);
pp = ones(size(sigma))*dist;
if edgefinger == 2
    pp = pp + sigma.*sin(phi);
end
end