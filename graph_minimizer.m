xl = [23 PG.S(end)];
yl = [-0.6 0.6];
% set(groot,'CurrentFigure',graphfig{2})
ax2 = graphfig{2}.Children;
set(ax2(1),'xlim',xl)
set(ax2(1),'ylim',yl)
set(ax2(3),'xlim',xl)
set(ax2(3),'ylim',yl)

%%
% xl = [12 PG.S(end)];
% yl = [-1 1.5];
% xl = [0 8];
% yl = [-0.4 2.1];
% xl = [0 32];
% yl = [0 2.5];
% xl = [44 52];
% yl = [-0.4 1.5];
xl = [4 8];
% yl = [-pi() pi()];
ax1 = graphfig{1}.Children;
set(ax1(1),'xlim',xl)
set(ax1(1),'ylim',yl)
set(ax1(3),'xlim',xl)
set(ax1(3),'ylim',yl)