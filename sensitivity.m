% parameters
r = 100;

b = 1000*r;
a = 0.5*b;
c = 0.5*b;
e = 2;
d = 5*e;
f = 2;
g = 0.001/r;
h = 0.0025/r;
k = 0.0015/r;
l = 0.01/r;

m=4.76;
T=100;

%u values
ulist=[0.000098,0.000099,0.0001,0.000101,0.000102];
%0.0000098,0.0000099,0.00001,0.0000101,0.0000102

for u=ulist

K1 = a*T^m;
K2 = b*T^m;
I0 = c*T^m;
r1 = d/T;
r2 = e/T;
r3 = f/T;
%alpha,beta,gamma,delta
n1 = g/T^(m+1);
n2 = h/T^(m+1);
n3 = k/T^(m+1);
n4 = l/T^(m+1);

%ODE function
func = @(t,y) [r1*y(1) - r1/K1*y(1)^2 - n1*y(1)*y(2)- n2*y(1)*y(3)
        r2*y(2) - r2/K2*y(2)^2 - n3*y(1)*y(2)
        -r3*y(3)-n4*y(1)*y(3)];
%lifespan
tf = T;
%initial condition
y0 = [u*K1 (1-u)*K2 I0];
%option-event function
options=odeset('Events',@event_func,'NormControl','on');
%run function
[t,y,te,ye,ie] = ode45(func, [0 tf], y0,options);

%graph
figure;
hold on;
plot(t,y(:,1));
plot(t,y(:,3));
hold off;

%display ratio
display(te/T);

end


%Event function to find intersection time
function [value,isterminal,direction]=event_func(t,y)
    value=y(1)-y(3);
    isterminal=0;
    direction=0;
end