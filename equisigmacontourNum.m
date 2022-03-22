function numsegments = equisigmacontourNum(PG,s1init,s2init)
syms t
[O1,O2] = PG.get('2Pos',s1init,s2init);
sigma = norm(O2-O1);
relvector = O2-O1;
Fullcircuit = 0;


s10 = s1init;
s20 = s2init;
theta0 = 0;
R = [cos(t) sin(t);-sin(t) cos(t)];
J = [0 -1;1 0];
thetavector = 0:-0.001:-pi();
thetadir = 1;
prevedge = [0;0];
startIn = 1;
path = [];
pathpoints = [];
s1vertex = [];
s2vertex = [];
s1vertexalt = [];
s2vertexalt = [];

while ~Fullcircuit %construct the contour from segments in each rectangle
% the contour passes through
try
    % Identify the edge-edge rectangle being explored
    noRotation = 0;
    edge1 = PG.get('edgeNum',s10);
    edge2 = PG.get('edgeNum',s20);
    if all([edge1;edge2]==prevedge) % check if the default rectangle is the one already previously visited by the loop
        %fix the direction of the contour by picking the neighboring
        %rectangle - setting the edge supported by the vert
        isvertex = any(PG.S(PG.VL)==s10)*1 + any(PG.S(PG.VL)==s20)*2;
        switch isvertex
            case 1
                edge1 = edge1-1;
                if edge1 == 0
                    edge1 = PG.nv;
                    s10 = PG.S(end);
                end
            case 2
                edge2 = edge2-1;
                if edge2 == 0
                    edge2 = PG.nv;
                    s20 = PG.S(end);
                end
            case 3 
                disp('vertex-vertex case')
                edge1 = edge1-1; %to prevent crossing a corner, we move s2 a tiny bit
                if edge1 == 0
                    edge1 = PG.nv;
                    s10 = PG.S(end);
                end
                if s20>0
                    s20 = s20-1E-6;
                else
                    s20 = s20+1E-6;
                end
        end
    end 
    if edge1 == 1 && s10 == PG.S(end)
        s10 = 0;
    end
    if edge2 == 1 && s20 == PG.S(end)
        s20 = 0;
    end   
    % find the functions describing the contour
    s1left = PG.S(PG.VL(edge1+1))-s10;
    s1leftalt = PG.S(PG.VL(edge1))-s10;
    s2left = PG.S(PG.VL(edge2+1))-s20;
    s2leftalt = PG.S(PG.VL(edge2))-s20;
    if any([s1left,s2left,s1leftalt,s2leftalt]==0)
        startIn = 0;
    end
    Rt0 = [cos(theta0) -sin(theta0);sin(theta0) cos(theta0)];
    n1 = Rt0*PG.normal(:,edge1);
    n2 = Rt0*PG.normal(:,edge2);
    if abs(det([J*n1 -J*n2]))>1E-10
        dsvector = [J*n1 -J*n2]\(R-eye(2))*relvector;
        Sv = thetavector;
        if abs(det([J*n1 -J*n2]))<0.5 %edges are close to parallel, small rotations cause high displacements
            % to have a high-resolution curve, add resolution to theta
            % accordingly
            Sv = thetavector(1):diff(thetavector(1:2))*10^max(-3,(round(log10(abs(det([J*n1 -J*n2]))),0))):thetavector(end);
        end
        % find the end of the segment in the current rectangle
        dsf = matlabFunction(dsvector);
        dsmatrix = feval(dsf,Sv);
        if startIn %special case when starting inside a rectangle
            s1vertex = find(dsmatrix(1,1:end)>=s1left,1,'first');
            s1vertexalt = find(dsmatrix(1,1:end)<=s1leftalt,1,'first');
            s2vertex = find(dsmatrix(2,1:end)>=s2left,1,'first');
            s2vertexalt = find(dsmatrix(2,1:end)<=s2leftalt,1,'first');
        else
            [s1vertex,s1vertexalt,s2vertex,s2vertexalt] = checkRectangle(dsmatrix,s1left,s1leftalt,s2left,s2leftalt,0);
            thetadir = 1;
            if all(isempty([s1vertex,s2vertex,s1vertexalt,s2vertexalt])) % need to rotate in the opposite direction
                Sv = -Sv;
                dsmatrix = feval(dsf,Sv);
                [s1vertex,s1vertexalt,s2vertex,s2vertexalt] = checkRectangle(dsmatrix,s1left,s1leftalt,s2left,s2leftalt,0);
                thetadir = -1;
            end
        end
        if isempty(s1vertex)
            s1vertex = inf;
        end
        if isempty(s2vertex)
            s2vertex = inf;
        end
        if isempty(s1vertexalt)
            s1vertexalt = inf;
        end
        if isempty(s2vertexalt)
            s2vertexalt = inf;
        end
        s1vertex = min([s1vertex,s1vertexalt]);
        s2vertex = min([s2vertex,s2vertexalt]);
        if min(s1vertex,s2vertex)==inf
            disp('intercept error');
        end
        Sv = Sv(1:min([s1vertex,s2vertex]));
        n1 = PG.normal(:,edge1);
        n2 = PG.normal(:,edge2);
        s1prev = s10;
        s2prev = s20;
        count = 0;
        while s1vertex==s2vertex %curve segments exits rectangle too close to corner to be clear which edge is crossed
            %add checking which edge is crossed
            Svtag = linspace(Sv(end-1),Sv(end),10);
            dsmatrixtag = feval(dsf,Svtag);
            [s1vertextag,s1vertexalttag,s2vertextag,s2vertexalttag] = checkRectangle(dsmatrixtag,s1left,s1leftalt,s2left,s2leftalt,1);
            if isempty(s1vertextag)
                s1vertextag = inf;
            end
            if isempty(s2vertextag)
                s2vertextag = inf;
            end
            if isempty(s1vertexalttag)
                s1vertexalttag = inf;
            end
            if isempty(s2vertexalttag)
                s2vertexalttag = inf;
            end
            mins1v = min([s1vertextag,s1vertexalttag]);
            mins2v = min([s2vertextag,s2vertexalttag]);
            minsv = min(mins1v,mins2v);
            s1vertex = s1vertex-2+mins1v;
            s2vertex = s2vertex-2+mins2v;
            svertex = min([s1vertex,s2vertex]);
            Sv = [Sv(1:end-2),Svtag(1:minsv)];
        end
        dsmatrix = feval(dsf,Sv);
        if s1vertex<s2vertex
            s10 = PG.S(PG.VL(edge1+(dsmatrix(1,min([s1vertex,s2vertex]))>0)));
            s21 = s20+dsmatrix(2,end);
            if s21>PG.S(end)
                s21 = s21-PG.S(end);
            end
            [O10,O20] = PG.get('2pos',s10,s21);
            sigmagal = norm(O20-O10);
            p = [n2 -J*n2]\(O10-O20);
            h = p(1);
            s20check = s21+sqrt(sigmagal^2-h^2)-sqrt(sigma^2-h^2);
            s20 = s21+sign(p(2))*(sqrt(sigmagal^2-h^2)-sqrt(sigma^2-h^2));
            O1check = O10;
            O2check = PG.get('1pos',s20);
            fingers= O2check-O1check;
            sigmacheck = norm(fingers)-sigma;
%             drawsigmaissue(PG,n2,sigma,s10,s21)
%             set(0, 'CurrentFigure', contspace);
            theta0new = atan2(relvector(2),relvector(1))-atan2(fingers(2),fingers(1));
            theta0prev = theta0;
            theta0 = theta0new;
        else
            s20 = PG.S(PG.VL(edge2+(dsmatrix(2,min([s1vertex,s2vertex]))>0)));
            s11 = s10+dsmatrix(1,end);
            if s11>PG.S(end)
                s11 = s11-PG.S(end);
            end
            [O10,O20] = PG.get('2pos',s11,s20);
            sigmagal = norm(O20-O10);
            p = [n1 -J*n1]\(O20-O10);
            h = p(1);
            s10 = s11+sign(p(2))*(sqrt(sigmagal^2-h^2)-sqrt(sigma^2-h^2));
            O1check = PG.get('1pos',s10);
            fingers= O20-O1check;
            sigmacheck = norm(fingers)-sigma;
%             drawsigmaissue(PG,n1,sigma,s20,s11)
%             set(0, 'CurrentFigure', contspace);
%             if abs(sigmacheck)>abs(sigmagal-sigma)
%                 disp('error start');
%             end
            theta0new = atan2(relvector(2),relvector(1))-atan2(fingers(2),fingers(1));
            theta0prev = theta0;
            theta0 = theta0new;
        end
        dsmatrix = [s1prev;s2prev]+dsmatrix;

    else % the supported edges are parallel
        
        dsvector = t*[1;sign(n1.'*n2)]; % the object isn't rotating, instead use ds1 for the parametrization of the segment
        if sign(n1.'*n2)>0
            Sv = [0:0.001:min([s1left,s2left]),min([s1left,s2left])];
        else
            Sv = [0:0.001:min([s1left,-s2leftalt]),min([s1left,-s2leftalt])];
        end
        if s1left==0
            if sign(n1.'*n2)>0
                Sv = [0:-0.001:max(s1leftalt,s2leftalt),max(s1leftalt,s2leftalt)];
            else
                Sv = [0:-0.001:max(s1leftalt,-s2left),max(s1leftalt,-s2left)];
            end
        end
        if s2left==0
            if sign(n1.'*n2)>0
                Sv = [0:-0.001:max(s1leftalt,s2leftalt),max(s1leftalt,s2leftalt)];
            else
                Sv = [0:0.001:min(s1left,-s2leftalt),min(s1left,-s2leftalt)];
            end
        end
        if s2leftalt==0
            if sign(n1.'*n2)>0
                Sv = [0:0.001:min(s1left,s2left),min(s1left,s2left)];
            else
                Sv = [0:-0.001:max(s1leftalt,-s2left),max(s1leftalt,-s2left)];
            end
        end
        if Sv(2)<0
            Sv = -Sv;
            dsvector = -dsvector;
        end
        dsf = matlabFunction(dsvector);
        dsmatrix = feval(dsf,Sv);
        noRotation = 1;
        s1prev = s10;
        s2prev = s20;
        s10 = s10+dsmatrix(1,end);
        s20 = s20+dsmatrix(2,end);
        dsmatrix = dsmatrix + [s1prev;s2prev];
        dsmatrix(:,end) = [s10;s20];
    end
    %     if theta0>2*pi()||theta0<-2*pi()
    %         Fullcircuit = 1;
    %     end
    if ~isempty(path)
        if find(all(path==[edge1;edge2],1))
            if find(all(abs(pathpoints-[s10;s20])<1E-4),1)
                Fullcircuit = 1;
            end
        end
    end
    path = [path,[edge1;edge2]];
    pathpoints = [pathpoints,[s10;s20]];
    
    
    if startIn % should only happen once, for the first segment
        startIn=0;
    end
    % add the resulting segment to the numerical and analytical
    % representations of the double-support contour
%     if ~noRotation
%         EquiSigmaNumeric = [EquiSigmaNumeric,[dsmatrix;theta0prev+Sv]];
%     else
%         EquiSigmaNumeric = [EquiSigmaNumeric,[dsmatrix;theta0*ones(size(Sv))]];
%     end
    
    % prepare the next iteration of the loop
    if imag(s10)~=0
        disp('big error')
    end
    hold on
    col = [0,sigma/20,1-sigma/20];
    if size(path,2)>60
        plot(dsmatrix(1,:),dsmatrix(2,:),'Color','r','linewidth',2)
%         disp('suspicious length')
    end
    plot(dsmatrix(1,:),dsmatrix(2,:),'Color',col,'linewidth',2)
    if ~any(abs(PG.S(PG.VL)-s10)>1E-14) && ~any(abs(PG.S(PG.VL)-s20)>1E-14)
        disp('something funky');
    end
%     if size(dsmatrix,2)<5
%          disp('weird length error');
%     end
    prevedge = [edge1;edge2];
    s1vertex = [];
    s2vertex = [];
    s1vertexalt = [];
    s2vertexalt = [];
catch errmsg
    errmsg.message
    errmsg.stack.line
    disp(errmsg)
end
end
numsegments = size(path,2);
end

function [s1vertex,s1vertexalt,s2vertex,s2vertexalt] = checkRectangle(dsmatrix,s1left,s1leftalt,s2left,s2leftalt,startIn)
s1vertex = [];
s2vertex = [];
s1vertexalt = [];
s2vertexalt = [];
% find on which of the rectangle's edge does the curve begin:
edge = 0; %edge number represents the edge in right-top-left-bottom order
if s1left ==0
    edge = edge+1;
end
if s2left ==0
    edge = edge+2;
end
if s1leftalt ==0
    edge = edge+3;
end
if s2leftalt ==0
    edge = edge+4;
end
% check if the curve found enters the right rectangle:
if ~startIn % startIn represents whether the curve segment starts inside the rectangle or on the boundary
    tangent = diff(dsmatrix(:,1:2),1,2); % if started on the boundary, check the direction of the curve to make sure it enters the right rectangle
    switch edge
        case 1
            if tangent(1)>0
                return
            end
        case 2
            if tangent(2)>0
                return
            end
        case 3
            if tangent(1)<0
                return
            end
        case 4
            if tangent(2)<0
                return
            end
    end
end
% check crossing the boundary only if it is not the starting boundary
% if s1left~=0
    s1vertex = find(dsmatrix(1,1:end)>s1left,1,'first');
%     if ~isempty(s1vertex)
%         ds2 = dsmatrix(2,s1vertex);
%         if ds2>s2left ||ds2<s2leftalt
%             s1vertex=[];
%         end
%     end
% end
% if s1leftalt~=0
    s1vertexalt = find(dsmatrix(1,1:end)<s1leftalt,1,'first');
%     if ~isempty(s1vertexalt)
%         ds2 = dsmatrix(2,s1vertexalt);
%         if ds2>s2left ||ds2<s2leftalt
%             s1vertexalt=[];
%         end
%     end
% end
% if s2left~=0
    s2vertex = find(dsmatrix(2,1:end)>s2left,1,'first');
%     if ~isempty(s2vertex)
%         ds1 = dsmatrix(1,s2vertex);
%         if ds1>s1left ||ds1<s1leftalt
%             s2vertex=[];
%         end
%     end
% end
% if s2leftalt~=0
    s2vertexalt = find(dsmatrix(2,1:end)<s2leftalt,1,'first');
%     if ~isempty(s2vertexalt)
%         ds1 = dsmatrix(1,s2vertexalt);
%         if ds1>s1left ||ds1<s1leftalt
%             s2vertexalt=[];
%         end
%     end
% end
end