% Author : Mathieu Bresciani
% More details on "Robust exact differentiation" : 
% http://www.tau.ac.il/~levant/mathappl/dfr_inst.pdf

clear all
close all
clc

Ts = 0.001;
t = 0:Ts:10;
nEpochs = length(t);
f = 0.5*sin(0.5*t)+0.5*cos(t);
f1 = 0.25*cos(0.5*t)-0.5*sin(t);

L = 1;

z0 = 0;
z1 = 0;

z0_mem = zeros(1,nEpochs);
z1_mem = zeros(1,nEpochs);
for k = 1:nEpochs
    d = z0-f(k);
    v0 = -1.5*L^(1/2)*abs(d)^(1/2)*sign(d)+z1;
    z0_dot = v0;
    z1_dot = -1.1*L*sign(d);
    
    z0 = z0 + z0_dot*Ts;
    z1 = z1 + z1_dot*Ts;
    
    z0_mem(k) = z0;
    z1_mem(k) = z1;
end

figure;
plot(t,f,t,z0_mem,t,f1,t,z1_mem);
legend('Position ref.','Position est.', 'True velocity', 'Velocity est')