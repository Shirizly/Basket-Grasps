halfInd = -1;
thetaSaddle = -pi()/3;
thetaNext = thetaSaddle+halfInd*pi();
thetaNextLimited = max(min(thetaNext,pi()),-pi());
thetaRange1 = sort([thetaSaddle, thetaNextLimited]);
thetaRange2 = [-2*pi(),-2*pi()];
if thetaSaddle*halfInd>0
ThetaEdge = -pi()*halfInd;
thetaRange2 = sort([ThetaEdge, thetaSaddle-halfInd*pi()]);
end
figure
hold on
plot(1:2,thetaRange1)
plot(2:3,thetaRange2)
ylim([-pi(),pi()])