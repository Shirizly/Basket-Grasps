function [Open_list,Closed_list,A,graph_nodes] = GraphSearch(PG,contv,graspind)
% 
    % PG.S,PG.X are discretization of the boundary
    % contv are the (s1,s2) of each point on the equisigma contour
    % graspind of the current basket grasp (in A)    
    
    % O and A are lists of neigbours, closed (checked) and open
    % NB is a matrix of neighborhoods
    plot_nodes = 0;
    %% define the positions of the fingers
    [basepos(:,1),basepos(:,2)] = getGraspPos(PG,contv(1,graspind),contv(2,graspind));
    %% find 2-finger max (saddle) and min nodes
    for i = 1:size(contv,2)
        [newpos(:,1),newpos(:,2)] = getGraspPos(PG,contv(1,i),contv(2,i));
        extpos(:,2*i-1:2*i) = newpos;
        [dsh(i),move] = PG.get('comRelativeHeight',basepos,newpos);
        theta(i) = move{1};
        d(:,i) = move{2};
    end
    smoothh = dsh; % smooth the graph since it is very noisy, in order to find peaks
    win = 10;
    for i = win+1:size(contv,2)-win
        smoothh(i) = mean(dsh(i-win:i+win));
    end
    for i = win+1:size(contv,2)-win
        smoothh(i) = mean(smoothh(i-win:i+win));
    end
    %% plot the height and find the peaks
%     figure
%     hold on
%     plot([graspind graspind],[min(h),max(h)],'--r')
%     plot(1:size(contv,2),smoothh);
    % plot(1:size(contv,2),h,'--c');
    
    [~,maxinds] = findpeaks(smoothh);
    [~,mininds] = findpeaks(-smoothh);
    removeind = [find(maxinds<win) find(maxinds>size(contv,2)-win)];
    maxinds(removeind) = [];
    removeind = [find(mininds<win) find(mininds>size(contv,2)-win)];
    mininds(removeind) = [];
%     plot(maxinds,smoothh(maxinds),'.r','markersize',10)
%     plot(mininds,smoothh(mininds),'.g','markersize',10)
%     grid on
%     hold off
%% find single support equilibrium: (the old way)
    % finding all the 2-finger grasps with a finger on the 1-finger grasp
    % point
%     [ss1,ssh]= PG.PG.SPG.SEqcheck(PG.S,PG.X,PG.VL);
%     ssequiv1 = []; ssequiv2 = []; ssequivh1 = []; ssequivh2 = [];
%     for i=1:size(ss1,2)
%         ssequiv1 = [ssequiv1 find(contv(1,:)==ss1(i))];
%         ssequivh1 = [ssequivh1 ones(length(find(contv(1,:)==ss1(i))))*(ssh(i)+basepos(2,1))];
%         ssequiv2 = [ssequiv2 find(contv(2,:)==ss1(i))];
%         ssequivh2 = [ssequivh2 ones(length(find(contv(2,:)==ss1(i))))*(ssh(i)+basepos(2,2))];
%     end
%     ssnodes = [ssequiv1 ssequiv2];
%     ssheight = [ssequivh1 ssequivh2];
    
%% find single support equilibrium: (the new way)
%     % finding 2-finger grasps that have the same com height as the
%     % 1-finger grasp
%     epsi = 2*sqrt(mean((diff(smoothh)).^2));
%     [sss,ssh]= PG.SSEqcheck(PG.S,PG.X,PG.VL);
%     ssequiv = []; ssequivh = [];
%     for j=1:size(sss,2)
%         h1 = ssh(j)+basepos(2,1);
%         ssind1 = find(abs(smoothh-h1)<epsi);
%         flag = length(ssind1)>0; i=1; rem = 0;
%         while flag % keeping only the real point out of a chain of close enough points
%             if ssind1(i+1)-ssind1(i)<=1+rem
%                 dif1 = abs(smoothh(ssind1(i))-h1);
%                 dif2 = abs(smoothh(ssind1(i+1))-h1);
%                 if dif1>dif2
%                     ssind1(i) = [];
%                     rem = 0;
%                 else
%                     ssind1(i+1) = [];
%                     rem = rem+1;
%                 end
%                 i=i-1;
%             else
%                 rem = 0;
%             end
%             i=i+1;
%             if i>=length(ssind1)-1
%                 flag  = false;
%             end
%         end
%         % check that the contact involves the relevant edge - function
%         ssind1 = checkedge(contv(1,:),PG.S(PG.VL),ssind1,sss(j));
%         
%         h2 = ssh(j)+basepos(2,2);
%         ssind2 = find(abs(smoothh-h2)<epsi);
%         flag = length(ssind2)>0; i=1; rem = 0;
%         while flag % keeping only the real point out of a chain of close enough points
%             if ssind2(i+1)-ssind2(i)<=1+rem
%                 dif1 = abs(smoothh(ssind2(i))-h1);
%                 dif2 = abs(smoothh(ssind2(i+1))-h1);
%                 if dif1>dif2
%                     ssind2(i) = [];
%                     rem = 0;
%                 else
%                     ssind2(i+1) = [];
%                     rem = rem+1;
%                 end
%                 i=i-1;
%             else
%                 rem = 0;
%             end
%             i=i+1;
%             if i>=length(ssind2)-1
%                 flag  = false;
%             end
%         end
%         % check that the contact involves the relevant edge - function
%         ssind2 = checkedge(contv(2,:),PG.S(PG.VL),ssind2,sss(j));
%         
%         ssind = [ssind1 ssind2];
%         ssequiv = [ssequiv ssind];
%         ssequivh = [ssequivh smoothh(ssind)];
%     end
%     ssnodes = ssequiv;
%     ssheight = ssequivh;
%% find single support equilibrium: (the new way)
    % finding 2-finger grasps that have the same com height as the
    % 1-finger grasp
    epsi = 2*sqrt(mean((diff(smoothh)).^2));
    [sss,ssh]= PG.SSEqcheck();
    ssequiv = []; ssequivh = [];
    for j=1:size(sss,2)
        h1 = ssh(j)+basepos(2,1);


%% Generate a struct of all graph nodes:
    indexlist = [maxinds,mininds,ssnodes];
    heightlist = [dsh(maxinds),dsh(mininds),ssheight];
    typelist = [ones(size(maxinds)),-1*ones(size(mininds)),zeros(size(ssnodes))];
    graph_nodes = struct('index',indexlist,'h',heightlist,'type',typelist);
    
    %% plot all nodes and their height
    if plot_nodes == 1
    hold on
    for i = 1:size(maxinds,2)
        plota(contv(:,maxinds(i)),'+r','markersize',10,'linewidth',3);
        textpos = contv(:,maxinds(i));
        textpos = textpos + (0.5*[2;0]*(textpos(1)<textpos(2))+1*[-2;0]*(textpos(1)>textpos(2)));
%         text(textpos(1),textpos(2),['h = ' num2str(round(dsh(maxinds(i)),2))],'fontweight','bold');
    end
    for i = 1:size(mininds,2)
        plota(contv(:,mininds(i)),'+g','markersize',10,'linewidth',3);
        textpos = contv(:,mininds(i));
        textpos = textpos + (0.5*[2;0]*(textpos(1)<textpos(2))+1*[-2;0]*(textpos(1)>textpos(2)));
%         text(textpos(1),textpos(2),['h = ' num2str(round(dsh(mininds(i)),2))],'fontweight','bold');
    end
    for i = 1:size(ssnodes,2)
        plota(contv(:,ssnodes(i)),'*r','markersize',10,'linewidth',3);
        textpos = contv(:,ssnodes(i));
        textpos = textpos + (0.5*[2;0]*(textpos(1)<textpos(2))+1*[-2;0]*(textpos(1)>textpos(2)));
%         text(textpos(1),textpos(2),['h = ' num2str(round(ssheight(i),2))],'fontweight','bold');
    end
    hold off
    end
    %%
    
    
    maxnb = size(graph_nodes.index,2);
    A = false(maxnb);
    flag = 1;
    mindiff = min(abs(graph_nodes.index-graspind));
    current_node = find(abs(graph_nodes.index-graspind)==mindiff);
    Open_list=current_node;
    Closed_list=[];
    while ~isempty(current_node)
        minnodes = find(graph_nodes.type==-1);
        edge_neighborup = minnodes(find(graph_nodes.index(minnodes)>graph_nodes.index(current_node),1,'first'));
        if isempty(edge_neighborup)
            firstmin = min(graph_nodes.index(minnodes));
            edge_neighborup = minnodes(find(graph_nodes.index(minnodes)==firstmin));
            neighborup = find(graph_nodes.index>graph_nodes.index(current_node));
            neighborup = [neighborup find(graph_nodes.index<graph_nodes.index(edge_neighborup))];
        else
            neighborup = find(and(graph_nodes.index>graph_nodes.index(current_node),graph_nodes.index<graph_nodes.index(edge_neighborup)));
        end
        edge_neighbordown = minnodes(find(graph_nodes.index(minnodes)<graph_nodes.index(current_node),1,'last'));
        if isempty(edge_neighbordown)
            lastmin = max(graph_nodes.index(minnodes));
            edge_neighbordown = minnodes(find(graph_nodes.index(minnodes)==lastmin));
            neighbordown = find(graph_nodes.index<graph_nodes.index(current_node));
            neighbordown = [neighbordown find(graph_nodes.index>graph_nodes.index(edge_neighbordown))];
        else
            neighbordown = find(and(graph_nodes.index<graph_nodes.index(current_node),graph_nodes.index>graph_nodes.index(edge_neighbordown)));
        end
        
        switch graph_nodes.type(current_node)
            case -1 % minimum nodes
                                neighbor = [neighborup,neighbordown];
            case {1,0}
                neighbor = [neighborup,edge_neighborup,neighbordown,edge_neighbordown];
        end
        
        A(current_node,neighbor) = true;
        % Add the neighbors of the current node to the open list:
        Open_list=union(Open_list,setdiff(find(A(current_node,:)).',[Open_list;Closed_list]));
        
        % Remove current node from open list:
        Open_list=setdiff(Open_list,current_node);
        
        % Add the current node to the end of the closed list:
        Closed_list=[Closed_list;current_node];
        
        % Find the next node with smallest sigma and define it as the new current
        % node:
        [~,min_ind]=min(graph_nodes.h(Open_list));
        current_node=Open_list(min_ind);
    end
    
end

