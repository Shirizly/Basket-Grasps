%% specific grasp contour - choosing a specific grasp
% 
[sf1,sf2] = ginput(1);
epsi = 1;
possiblex = find(abs(s1-sf1)<epsi);
possibley = find(abs(s2-sf2)<epsi);
possiblegrasp = intersect(possiblex,possibley);
mindist = 2*epsi^2;
for i=1:length(possiblegrasp)
    dist = (s1(possiblegrasp(i))-sf1)^2+(s2(possiblegrasp(i))-sf2)^2;
    if mindist>dist
        mindist=dist;
        grasp = possiblegrasp(i);
    end
end
sf1 = s1(grasp);
sf2 = s2(grasp);
hold on % adding a marking of the grasp to the contact-space plot
plot(sf1,sf2,'.k','MarkerSize',12);
hold off
[f1pos,f2pos] = PG.get('2Pos',sf1,sf2);
%% marking the specific grasp on the polygon
fig = PG.drawPolyGrasp(sf1,sf2,f1pos,f2pos);

%% Identifying the grasp's distance and contour
[cont,graspind,sig] = PG.GetSigmaContour(Sigma,sf1,sf2);
%% Plot the resulting contour
figure
contour(S,S,Sigma,c_layers/2)
xlabel('s_1')
ylabel('s_2')
colorbar
hold on
for i =1:numel(cont)
plot(cont{i}(1,:),cont{i}(2,:),'k','linewidth',1)
end
plot(sf1,sf2,'.k','MarkerSize',14)
% plot(cont{1}(1,graspind),cont{1}(2,graspind),'+k','MarkerSize',12)
% axis([0,maxS,0,maxS])
base = [sf1,sf2]; % original grasp
[basepos(:,1),basepos(:,2)] = PG.get('2Pos',base(1),base(2));
contv = cell2mat(cont.');
