%% ========================================================================
% ---- ACS6101 Assignment week 4
% ---- Name: Paulo Roberto Loma Marconi
% ---- 23/10/2018
%% === Question 1 =========================================================
clear; clc; close all;
%% plant G parameters
s = tf('s');
a = 2.5; b = 27;
%% a) Gain K for PO=10
PO = 10; % percentage overshoot
zeta = log(100/PO)/sqrt(pi^2+ (log(100/PO))^2 ); % damping ratio
PM_d = round(100*zeta)+1; % PM desired at the nearest round value
h = tand(180-90-PM_d); 
omega_c = roots([1 (a+b)/h -a*b]); % new omega_c
K = omega_c(2)*sqrt(omega_c(2)^2+a^2)*sqrt(omega_c(2)^2+b^2); % new gain K
G = K/( s*(s+a)*(s+b) ); % Plant G
figure(1);
margin(G);
figure(2);
step(feedback(G,1));
%% b) Phase-lead compensator
Kv = 25; % velocity error
Kg = Kv*a*b; % loop gain to satisfy Kv
G = Kg/( s*(s+a)*(s+b) ); % Plant G
% omega_c for the uncompensated system
omega_c = real( sqrt( roots([1 a^2+b^2 a^2*b^2 -Kg^2]) ) );
% obtaining actual PM
PM_act = 180-90-atand(omega_c(3)/a)-atand(omega_c(3)/b);
% additional phase angle from the compensator PM
theta = 0; % factor of correction
PM_c = round(PM_d+theta-PM_act);
% calculating alpha
alpha = round( (sind(PM_c)+1)/(1-sind(PM_c)) );
% Determine the new cross over frequency
omega_c_new = real( sqrt( roots([1 a^2+b^2 a^2*b^2 -Kg^2*alpha]) ) );
% Calculating the zero of the compensator
z = omega_c_new(3)/sqrt(alpha);
% Calculating the pole of the compensator
p = alpha*z;
% lead compensator
Gc = (s+z)/(s+p);
% new gain of the compensated system
K_new = sqrt(alpha)*Kg;
% compensated open-loop 
Gol = K_new*Gc*G/Kg;
% closed-loop 
Gcl = feedback(Gol,1);
%% c) evaluating the system to unit step and unit ramp
figure(2);
margin(Gol);
figure(3);
subplot(2,1,1);
step = 1/s; % step input
impulse(step,Gcl*step); 
title('Step response');
subplot(2,1,2);
ramp = 1/s^2; % ramp input
impulse(ramp,Gcl*ramp);
title('Ramp response');
% evaluating Kv for the new compensated system
Kv_new = (z/p)*K_new/(a*b);
% steady-state error to a unity ramp
ess_ramp = 1/Kv_new;
stepinfo(Gcl)
% bandwidth of Gcl
BW = bandwidth(Gcl);
% resonant frequency
omega_r = 4/(0.6*1.03) *sqrt(1-2*0.6^2);

% syms s
% ess = double(limit(s*Gol,s,0))