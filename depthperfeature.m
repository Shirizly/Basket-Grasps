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
if length(sigmaseg)==1
    N = 1;
end
sigma = linspace(minsig,maxsig,N); % segment of relevant sigmas
phi = phiseg;
if N>1
phi = interpn(sigmaseg,phiseg,sigma);
end
if vertfinger == 2 %switch the angles to work around finger 2, fix in the end
    phi = phi+pi();
end
% find puncture point com position
if vertnum>PG.nv
    vertnum = 1;
end
if edgenum>PG.nv
    edgenum = 1;
end
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
edge2 = featurepair(2);
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
        if length(sigmaseg)==1
            N = 1;
        end
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