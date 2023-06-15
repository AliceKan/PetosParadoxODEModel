%parameters
r = 10;

b = 1000*r;
a = 0.5*b;
c = 0.5*b;
e = 2;
d = 5*e;
f = 2;
g = 0.001/r;
h = 0.0015/r;
k = 0.0012/r;
l = 0.01/r;

figure;
T = 1;

m=4.76;

K1 = a*T^m;
K2 = b*T^m;
I0 = c*T^m;
r1 = d/T;
r2 = e/T;
r3 = f/T;
n1 = g/T^(m+1);
n2 = h/T^(m+1);
n3 = k/T^(m+1);
n4 = l/T^(m+1);

%ode system
func = @(t,y) [r1*y(1) - r1/K1*y(1)^2 - n1*y(1)*y(2)- n2*y(1)*y(3)
        r2*y(2) - r2/K2*y(2)^2 - n3*y(1)*y(2)
        -r3*y(3)-n4*y(1)*y(3)];


%phase portrait
subplot(3,1,[1,2]);
x=linspace(0,K1,5);
y=linspace(0,K2,10);
z=linspace(0,I0,5);

[X,Y,Z]=meshgrid(x,y,z);
u=zeros(size(X));
v=zeros(size(Y));
w=zeros(size(Z));


for i=1:numel(u)
    
    Yprime=func(0,[X(i);Y(i);Z(i)]);
    u(i)=Yprime(1);
    v(i)=Yprime(2);
    w(i)=Yprime(3);
    
    Vmod=sqrt(u(i)^2+v(i)^2+w(i)^2);
    u(i)=u(i)/Vmod;
    v(i)=v(i)/Vmod;
    w(i)=w(i)/Vmod;
end

quiver3(X,Y,Z,u,v,w);

xlabel('C');
ylabel('H');
zlabel('I');

hold on;
tf = T;
y0 = [0.0001*K1 0.9999*K2 I0];
opts=odeset('NormControl','on');
[t,y] = ode45(func, [0 tf], y0);
%color with time
color = linspace(0,1,length(t));
steps = t;

scatter3(y(:,1), y(:,2), y(:,3),30, steps, 'filled','MarkerFaceAlpha', 0.7);
colormap(jet(length(t)));
colorbar;
hold off;

%line graph
subplot(3,1,3);
hold on;
plot(t,y(:,1));
plot(t,y(:,2));
plot(t,y(:,3));
hold off;
legend('C','H','I');
xlabel('time');
ylabel('cell counts');




