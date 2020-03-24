function [ sigma_matrix ] = inter_finger_distance( X1,X2)
% inter_finger_distance computes a matrix containing the distance between 
% point pairs on the boundary of the polygonal object, corresponding to the 
% contact space coordinates, sigma=||X(s1)-X(s2)||.
%%
sigma_matrix=zeros(length(X1(:,1)));

for j=length(X1(:,1)):-1:1
    for i=1:length(X1(:,1))
        sigma_matrix(j,i)=sqrt((X2(j,1)-X1(i,1))^2+(X2(j,2)-X1(i,2))^2);
    end
end