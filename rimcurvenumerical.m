function rim = rimcurvenumerical(PG,ContactHeight,U,PP)
theta = -pi():0.01:pi();
[ height_matrix ] = single_support_height(PG,PG.S,theta);
height_matrix = height_matrix + ContactHeight;
c_layers=[U,U];
cont = contour(PG.S,theta,height_matrix.',c_layers);

rim = cont
end