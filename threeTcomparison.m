%parameters
r = 100;

b = 1000*r;
a = 0.5*b;
c = 0.5*b;
e = 1000;
d = 5*e;
f = 2;
g = 4.6/r;
h = 1.6/r;
k = 1.5/r;
l = 0.01/r;

figure;
T = 1;
i = 1;
while i <= 3

    K1 = a*T^4.76;
    K2 = b*T^4.76;
    I0 = c*T^4.76;
    r1 = d/T;
    r2 = e/T;
    r3 = f/T;
    n1 = g/T^5.76;
    n2 = h/T^5.76;
    n3 = k/T^5.76;
    n4 = l/T^5.76;

    %----omit the following for original ode equation----
    r10=r1*T;
    r20=r2*T;
    r30=r3*T;
    n10=n1*K2*T;
    n20=n2*I0*T;
    n30=n3*K1*T;
    n40=n4*K1*T;
    %----------------------------------------------------

    %nondimensionalized equations system
    func = @(t,y) [r10*y(1)*(1-y(1)) - n10*y(1)*y(2)- n20*y(1)*y(3)
        r20*y(2)*(1-y(2)) - n30*y(1)*y(2)
        -r30*y(3)-n40*y(1)*y(3)];
    tf = 1;
    y0 = [0.0001 0.9999 1];



    %---use the following instead for original ode equation----
    %func = @(t,y) [r1*y(1)*(1-y(1)/K1) - n1*y(1)*y(2)- n2*y(1)*y(3)
    %    r2*y(2)*(1-y(2)/K2) - n3*y(1)*y(2)
    %    -r3*y(3)-n4*y(1)*y(3)];
   
    %tf = T;
    %y0 = [0.0001*K1 0.9999*K2 I0];
    %----------------------------------------------------------

    opts=odeset('NormControl','on');

    [t,y] = ode45(func, [0 tf], y0,opts);

    subplot(3,1,i);
    plot(t,y);
    legend('C','H','I');
    xlabel('nondimensionalized time');
    ylabel({'nondimensionalized'; 'cell counts'});

    %also please change labels accordingly if using original ode equations

    T = T*10;
    i = i+ 1;
    
end