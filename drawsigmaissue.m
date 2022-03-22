function drawsigmaissue(PG,nu,sigma,sk,su)
J = [0 -1;1 0];
[Ok,Ou] = PG.get('2pos',sk,su);
sigmagal = norm(Ok-Ou);
p = [nu -J*nu]\(Ok-Ou);
h = p(1);
sm = su+sign(p(2))*(sqrt(sigmagal^2-h^2)-sqrt(sigma^2-h^2));
Om = PG.get('1pos',sm);
fingers= Ok-Om;
sigmacheck = norm(fingers)-sigma;

figure
hold on
plot([Ou(1),Om(1)],[Ou(2),Om(2)],'k')
plot(Om(1),Om(2),'r+')
plot(Ou(1),Ou(2),'g+')
plot([Ou(1),Ou(1)+nu(2)*p(2)],[Ou(2),Ou(2)-nu(1)*p(2)],'y--')
plot([Ok(1),Ok(1)+h*(-nu(1))],[Ok(2),Ok(2)+h*(-nu(2))],'g')
theta = 0:0.00001:2*pi();
plot(Ok(1)+sigma*cos(theta),Ok(2)+sigma*sin(theta),'b')
plot(Ok(1),Ok(2),'r+')
axis equal
hold off
end