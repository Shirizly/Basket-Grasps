function varargout = get(PG,spec,varargin)
% this function returns a parameter (or multiple) of the polygon, according
% to the spec, and according to additional information held in vargin
switch spec
    case 'edge'
        edgeNum = varargin{1};
        v1 = PG.vertex(:,edgeNum);
        if edgeNum==PG.nv
            v2 = PG.vertex(:,1);
        else
            v2 = PG.vertex(:,edgeNum+1);
        end
        normal = PG.normal(:,edgeNum);
        varargout = {v1,v2,normal};
        
        
    case 'normal'
        s = varargin{1};
        for k=1:PG.nv
            if s <= PG.sv(k+1)
                varargout{1} = PG.normal(:,k);
                return
            end
        end
        varargout{1} = [];
        
        
    case 'rect'
        sf1 = varargin{1};
        sf2 = varargin{2};
        s1min = PG.S(PG.VL(find(PG.S(PG.VL)<sf1,1,'last')));
        s2min = PG.S(PG.VL(find(PG.S(PG.VL)<sf2,1,'last')));
        s1max = PG.S(PG.VL(find(PG.S(PG.VL)>sf1,1,'first')));
        s2max = PG.S(PG.VL(find(PG.S(PG.VL)>sf2,1,'first')));
        varargout{1} = [s1min s1max s2min s2max];
        
        
    case '2Pos' %interpolate the original location of 2 points on the boundary
        s1 = varargin{1};
        s2 = varargin{2};
        ind12 = find(PG.S>=s1,1,'first');
        ind11 = find(PG.S<=s1,1,'last');
        if ind12~=ind11
            alpha1 = (PG.S(ind12)-s1)/(PG.S(ind12)-PG.S(ind11));
            f1pos = (alpha1*PG.X(ind11,:)+(1-alpha1)*PG.X(ind12,:)).';
        else
            f1pos = PG.X(ind12,:).';
        end % in case the points are the same, no need (or possible) to interpolate
        ind22 = find(PG.S>=s2,1,'first');
        ind21 = find(PG.S<=s2,1,'last');
        if ind22~=ind21
            alpha2 = (PG.S(ind22)-s2)/(PG.S(ind22)-PG.S(ind21));
            f2pos = (alpha2*PG.X(ind21,:)+(1-alpha2)*PG.X(ind22,:)).';
        else 
            f2pos = PG.X(ind22,:).';
        end
        varargout = {f1pos,f2pos};
        
        
    case '1Pos' %interpolate the original location of a point on the boundary
        s = varargin{1};
        ind2 = find(PG.S>=s,1,'first');
        ind1 = find(PG.S<=s,1,'last');
        if ind2~=ind1
        alpha = (PG.S(ind2)-s)/(PG.S(ind2)-PG.S(ind1));
        fpos = (alpha*PG.X(ind1,:)+(1-alpha)*PG.X(ind2,:)).';
        else % in case the points are the same, no need (or possible) to interpolate
            fpos = PG.X(ind2,:).';
        end
        varargout = {fpos};
        
        
    case 'edgeNum'
        s = varargin{1};
        edgenum = find(PG.S(PG.VL)>s,1,'first') - 1;
        if isempty(edgenum) %if point is last s (representing first vertex)
            edgenum = 1;
        end
        varargout{1} = edgenum;
        
    case 'comRelativeHeight'
        basepos = varargin{1};
        newpos = varargin{2};
        fingerline =basepos(:,2)-basepos(:,1);
        fl = fingerline/norm(fingerline);
        newfl = newpos(:,2) - newpos(:,1);
        newfl = newfl/norm(newfl);
        cs = [newfl(1) newfl(2);newfl(2) -newfl(1)]\fl;
        theta = atan2(cs(2),cs(1)); % angle clockwise
        d = basepos(:,1) - newpos(:,1);        
        R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        CoM = newpos(:,1) + R*(PG.com-newpos(:,1))+d;
        dist = CoM-PG.com;
        h = dist(2,1);
        varargout{1} = h;
        switch nargout
            case 2
                varargout{2} = {theta,d};
            case 3
                varargout{2} = R;
                varargout{3} = d;
        end
        
        
        case 'comLocationDS' % receives basepos as the current grasp, s1,s2 as the locations on the boundary
            % returns the location of the com, 
        basepos = varargin{1};
        s1 = varargin{2};
        s2 = varargin{3};
        [fs1,fs2] = PG.get('2Pos',s1,s2);
        newpos = [fs1,fs2];
        fingerline =basepos(:,2)-basepos(:,1);
        fl = fingerline/norm(fingerline);
        newfl = newpos(:,2) - newpos(:,1);
        newfl = newfl/norm(newfl);
        cs = [newfl(1) newfl(2);newfl(2) -newfl(1)]\fl;
        theta = -atan2(cs(2),cs(1)); % angle counterclockwise
        d = basepos(:,1) - newpos(:,1);        
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        CoM = newpos(:,1) + R*(PG.com-newpos(:,1))+d;
        dist = CoM;
        h = dist(2,1);
        varargout{1} = h;
        switch nargout
            case 2
                varargout{2} = {theta,d};
            case 3
                varargout{2} = R;
                varargout{3} = d;
        end
        
        
    case 'comLocationSS'
        f1 = varargin{1};
        s1 = varargin{2};
        theta = varargin{3};
        relative_to_contact = [cos(theta) -sin(theta);sin(theta) cos(theta)]*(PG.com-PG.get('1Pos',s1).');
        varargout{1} = f1+relative_to_contact;
end
        
end

%% old functions
% function [v1,v2,normal] = getEdge(PG,edgeNum)
%         v1 = PG.vertex(:,edgeNum);
%         if edgeNum==PG.nv
%             v2 = PG.vertex(:,1);
%         else
%             v2 = PG.vertex(:,edgeNum+1);
%         end
%         normal = PG.normal(:,edgeNum);
%     end
%     
%     function nm = getnormal(PG,s)
%         nm = [];
%         for k=1:PG.nv
%             if s <= PG.sv(k+1)
%                 nm = PG.normal(:,k);
%                 return
%             end
%         end
%     end
%     
%     function rect = getRectangle(PG,sf1,sf2)
%         s1min = PG.S(PG.VL(find(PG.S(PG.VL)<sf1,1,'last')));
%         s2min = PG.S(PG.VL(find(PG.S(PG.VL)<sf2,1,'last')));
%         s1max = PG.S(PG.VL(find(PG.S(PG.VL)>sf1,1,'first')));
%         s2max = PG.S(PG.VL(find(PG.S(PG.VL)>sf2,1,'first')));
%         rect = [s1min s1max s2min s2max];
%     end
%     
%     function [f1pos,f2pos] = getGraspPos(PG,s1,s2)
% %         ind1 = find(PG.S==s1);
%         ind1 = find(PG.S>=s1,1,'first');
% %         ind2 = find(PG.S==s2);
%         ind2 = find(PG.S>=s2,1,'first');
%         f1pos = PG.X(ind1,:);
%         f2pos = PG.X(ind2,:);        
%     end
%     
%     function f1pos = getSPos(PG,s1)
%         ind1 = find(PG.S>=s1,1,'first');
%         f1pos = PG.X(ind1,:);      
%     end
%     function edge_num = getEdgeNum(PG,s1)
%         edge_num = find(PG.S(PG.VL)>s1,1,'first') - 1;
%     end
%     
%     function [h,varargout] = getCMHeight(PG,basepos,newpos)
%         fingerline =basepos(:,2)-basepos(:,1);
%         fl = fingerline/norm(fingerline);
%         newfl = newpos(:,2) - newpos(:,1);
%         newfl = newfl/norm(newfl);
%         cs = [newfl(1) newfl(2);newfl(2) -newfl(1)]\fl;
%         theta = atan2(cs(2),cs(1)); % angle clockwise
%         % R = [cos(theta) sin(theta); -sin(theta) cos(theta)]; % self check
%         % R*newfl
%         d = basepos(:,1) - newpos(:,1);        
%         R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
%         CoM = newpos(:,1) + R*(PG.com-newpos(:,1))+d;
%         dist = CoM-PG.com;
%         h = dist(2,1);
%         switch nargout
%             case 2
%                 varargout{1} = {theta,d};
%             case 3
%                 varargout{1} = R;
%                 varargout{2} = d;
%         end
%             
%     end
%     
%     function com_pos = getCOMLoc(PG,f1,s1,theta) %find location of COM in a ss grasp
%         % body is held in a single-support grasp, its angle is theta
%         % (relative to original), finger's position is f1, contact point on
%         % the boundary is s1, output is location of CoM and it's rotation;
%         relative_to_contact = [cos(theta) -sin(theta);sin(theta) cos(theta)]*(PG.com-PG.getSPos(s1).');
%         com_pos = f1+relative_to_contact;
%     end
