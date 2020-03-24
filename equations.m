h = 1.5;
L = 1;
dir = -1;
a = dir*(0:0.001:L);


t = (pi() - 2*atan2(h,-a));

y = h*cos(t)-a.*sin(t);
x = -a.*cos(t)-h*sin(t);

hold on
plot(t,x,'k')

%% % % --------------------

tstart = t(end);
% h = 0.1;
x0 = x(end);
L = 1;
t0 = pi()/2+tstart;
phi0 = atan2(h,x0);
df = phi0-t0;

% syms a
r = (h^2+x0^2)^0.5;

b = r*cos(df);
disc = b^2-x0^2;
alpend1 = b-sqrt(disc);
% alpend2 = b+sqrt(disc);
if disc<0
    alpend1 = L;
end
alpend = min(alpend1,L);

a = dir*[0:0.0001:alpend alpend];


c = (r^2+a.^2-2*r*a*cos(df)).^0.5;
b = atan2(r*sin(df),r*cos(df)-a);
% beta = atan2(-a*sin(-df),-a*cos(-df)+r);
tgal = asin(h./c)-b;
% phi = asin(h./c)-beta;
y = r*sin(tgal+df)-a.*sin(tgal);
x = -a.*cos(tgal)+r*cos(tgal+df);

t = tgal-(t0-tstart);
% rem = find(imag(t)~=0,1,'first');
% 
% t(rem) = [];
% a(rem) = [];
% x(rem) = [];
% y(rem) = [];
hold on
plot(t,x,'k')

%%
for i=1:-1
    
% h = 2.5;
% L = 2;
a = dir*(0:0.001:L);


t = zeros(size(a));

y = a.*sin(t);
x = -a;

hold on
plot(t,x,'k')

% % % - - - - - - - - - - -

tstart = t(end);
% h = 0.1;
x0 = x(end);
L = 1;
t0 = pi()/2+tstart;
phi0 = atan2(h,x0);
df = phi0-t0;

% syms a
r = (h^2+x0^2)^0.5;

b = r*cos(df);
disc = b^2-x0^2;
alpend1 = b-sqrt(disc);

if disc<0
    alpend1 = L;
end
alpend = min(alpend1,L);

a = [0:0.0001:alpend alpend];

c = (r^2+a.^2-2*r*a*cos(df)).^0.5;
beta = atan2(r*sin(df),r*cos(df)-a);

tgal = pi()-asin(h./c)-beta;

y = r*sin(tgal+df)-a.*sin(tgal);
x = -a.*cos(tgal)+r*cos(tgal+df);

t = tgal-(t0-tstart);

hold on
plot(t,x,'k')
end
%%
tstart = t(end);
% h = 0.1;
x0 = x(end);
L = 1;
t0 = +tstart;
phi0 = atan2(h,x0);
df = phi0-t0;

% syms a
r = (h^2+x0^2)^0.5;

a = -a;


c = (r^2+a.^2-2*r*a*cos(df)).^0.5;
beta = atan2(r*sin(df),r*cos(df)-a);

tgal = pi()-asin(h./c)-beta;

y = r*sin(tgal+df)-a.*sin(tgal);
x = -a.*cos(tgal)+r*cos(tgal+df);

t = tgal-(t0-tstart);

hold on
plot(t,x,'k')




