function [BGDepth,BGSigma,EscapeFeatures] = DSEscape(PG,BGSeg,BGSigma,BGPhi)
%run the computation once for each direction
[BGDepth1,~,EscapeFeatures1] = DSEscape1(PG,BGSeg,BGSigma,BGPhi,1,1,cell(1,2));
[BGDepth2,~,EscapeFeatures2] = DSEscape1(PG,BGSeg,BGSigma,BGPhi,1,-1,cell(1,2));
% afterwards, combine to take the minimum of the two double-support graphs
[BGDepth,Index] = min([BGDepth1;BGDepth2],[],1);
EscapeFeatures = EscapeFeatures1;
EscapeFeatures(1,:) = EscapeFeatures1(1,:).*(Index==1) + EscapeFeatures2(1,:).*(Index==2);
EscapeFeatures(2,:) = EscapeFeatures1(2,:).*(Index==1) + EscapeFeatures2(2,:).*(Index==2);
end
function [BGDepth,BGSigma,EscapeFeatures] = DSEscape1(PG,BGSeg,BGSigma,BGPhi,steps,direction,DataSet)
[bg,~,~] = bgheight(PG,BGSeg(:,[1,end]));
CurrentSRange = BGSeg(:,[1,end]);
BGDepth = inf*ones(size(BGSigma));
EscapeFeatures = zeros(2,size(BGSigma,2));
steps = steps-1;
while all(BGDepth == inf)
    steps = steps+1;
%      tic
    [EquiSigContours,DataSet] = equisigmacontourSym2(PG,CurrentSRange,steps,direction,DataSet);
%     toc
    possiblefeaturelist = zeros(2,2*numel(EquiSigContours));
    for i = 1:numel(EquiSigContours)
        % find which features are "visited" by the equiSigmaContours
        % each element is a sub-segment with the same feature pair
        EquiSigContour = EquiSigContours{i};
        EquiSigContourEndpoint = EquiSigContour{end-2};
        EquiSigSegment = EquiSigContourEndpoint{end};
        edgeedgepair = EquiSigSegment{2}.';
        vertexedgepair = EquiSigSegment{end}.';
        possiblefeaturelist(:,2*i-1:2*i) = [edgeedgepair,vertexedgepair];
    end
    possiblefeaturelist = unique(possiblefeaturelist.','rows').';
%     timechar = toc
%     tic
    for i=1:size(possiblefeaturelist,2)
        featurepair = possiblefeaturelist(:,i);
        %         try
        [ppv,sigma] = depthperfeature2(PG,BGSigma,BGPhi,featurepair);
        %         catch ME
        %             disp('what?')
        %         end
        interpbg = linspace(bg(2,1),bg(2,2),length(sigma));
        [BGDepth,Index] = min([BGDepth;ppv-interpbg],[],1);
        EscapeFeatures(:,Index==2) = [featurepair(1)*ones(1,nnz(Index==2));featurepair(2)*ones(1,nnz(Index==2))];
%         disp(BGSigma([find(Index==2,1,'first'),find(Index==2,1,'last')]))
        if any(BGDepth<-1E-3)
            disp('theoretical error!');
        end
    end
%     timecomp = toc
end
complete = ~any(BGDepth==inf);

while ~complete
    cutoffpoint1 = find(BGDepth==inf,1,'first');
    cutoffpoint2 = cutoffpoint1 + find(BGDepth(cutoffpoint1+1:end)~=inf,1,'first')-1;
    if isempty(cutoffpoint2)
        cutoffpoint2 = length(BGDepth);
    end
    BGSubSeg = BGSeg(:,cutoffpoint1:cutoffpoint2);
    BGSubSigma = BGSigma(cutoffpoint1:cutoffpoint2);
    BGSubPhi = BGPhi(cutoffpoint1:cutoffpoint2);
    [BGSubDepth,~,featurepairSub] = DSEscape1(PG,BGSubSeg,BGSubSigma,BGSubPhi,steps+1,direction,DataSet);
    BGDepth(cutoffpoint1:cutoffpoint2) = min([BGDepth(cutoffpoint1:cutoffpoint2);BGSubDepth],[],1);
    EscapeFeatures(:,cutoffpoint1:cutoffpoint2) = featurepairSub;
    complete = ~any(BGDepth==inf);
end

end

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

% function BGSegPoint = BGFromSigma(PG,BGSegEnds,BGSegEndsSigma,BGSigmaPoint)
% %receives as input the current si endpoints and their sigma value, as well
% %as a sigma value inside the range.
% %returns the si values related to the sigma value given
% BGSEndpoint1 = BGSegEnds(:,1);
% BGSEndpoint2 = BGSegEnds(:,2);
% alpha = (BGSigmaPoint-BGSegEndsSigma(1))/(BGSegEndsSigma(2)-BGSegEndsSigma(1));
% BGSegPoint = BGSEndpoint1+alpha*(BGSEndpoint2-BGSEndpoint1);
% BGP1 = PG.get('1pos',BGSegPoint(1));
% BGP2 = PG.get('1pos',BGSegPoint(2));
% BGPrel = BGP2-BGP1;
% sigcheck = norm(BGPrel);
% if abs(sigcheck-BGSigmaPoint)>1E-5
%     disp('calc sigma error');
% end
% end
