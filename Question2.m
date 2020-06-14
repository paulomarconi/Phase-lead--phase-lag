%% ========================================================================
% ---- ACS6101 Assignment week 4
% ---- Name: Paulo Roberto Loma Marconi
% ---- 27/10/2018
%% === Question 2 =========================================================
clear; clc; close all;
%% plant G parameters
s = tf('s');
a = 9; b = 50;
%% a) Location of the dominant poles for PO=20 and ts=2.9
PO = 20; % percentage overshoot
ts = 2.9; % settling time
zeta = log(100/PO)/sqrt(pi^2+ (log(100/PO))^2 ); % damping ratio
zeta = round(zeta,2)-0.01;
omega_n = 4/(zeta*ts);
% desired location of dominant poles
s1 = -omega_n*zeta+omega_n*sqrt(1-zeta^2)*1i;
%% c) Phase lead compensator
z_lead = 1; % location of the desired zero
% using the angle condition to determine the location of the pole
x = -real(s1); y = imag(s1);
h = 180 + 180 - atand(y/(x-z_lead))-2*( 180-atand(y/x) )-atand(y/(a-x))...
    -atand(y/(b-x));
p_lead = y/tand(h)+x;
% Gc compensator TF
Gc_lead = (s+z_lead)/(s+p_lead);
% obtaining the gain K with the gain condition equation
K = ((-x)^2+y^2)*sqrt((-x+a)^2+y^2)*sqrt((-x+b)^2+y^2)*...
    sqrt((-x+p_lead)^2+y^2)/sqrt((-x+z_lead)^2+y^2);
% the open-loop system
G = K/(s^2*(s+a)*(s+b));% plant G
Gol_lead = Gc_lead*G; % open-loop with the Lead compensator
Gcl_lead = feedback(Gol_lead,1); % closed-loop with the Lead compensator

fig = figure(1);
rlocus(Gcl_lead);
hold;
% ploting the s1 and zeta in the rlocus
n = 0:1:160; m = n*sqrt(zeta^2/(1-zeta^2));
axis ([ -4 1 -4 4]);
plot (-m,n,'--'); % zeta
plot (-x,y,'rd');
saveas(fig,'Q2_Lead_rlocus.png');

fig = figure(2);
step(Gcl_lead);
saveas(fig,'Q2_Lead_step.png');

fig = figure(3);
margin(Gcl_lead);
saveas(fig,'Q2_Lead_margin.png');
BW_lead = bandwidth(Gcl_lead); % bandwidth

%% Using prefilter to reduce the overshoot
pf = z_lead; % selecting the zero of the lead compensator (z_lead)
Gpf = pf/(s+pf);

fig = figure(3);
step(Gpf*Gcl_lead);
saveas(fig,'Q2_Lead+prefilter_step.png');

fig = figure(4);
margin(Gpf*Gcl_lead);
saveas(fig,'Q2_Lead+prefilter_margin.png');
BW_lead_prefilter = bandwidth(Gpf*Gcl_lead); % bandwidth

%% d) Phase-lag compensator
ess_a = 0.025; % steady-state error for a parabolic input 0.5At^2
Ka_d = 1/ess_a; % Ka desired
Ka_act = (z_lead/p_lead)*(K/(a*b));
% calculating the zero and pole of the compensator
alpha = Ka_d/Ka_act;
z_lag = abs(x)/10;
p_lag = z_lag/alpha;
% Gc phase-lag compensator
Gc_lag = (s+z_lag)/(s+p_lag);
%% Evaluating the performance with the phase-lead, phase-lag and prefilter
Gol_lead_lag = Gc_lead*Gc_lag*G;
Gcl_lead_lag = feedback(Gol_lead_lag,1);

fig = figure(5); 
rlocus(Gcl_lead_lag);
hold;
% ploting the s1 and zeta in the rlocus
n = 0:1:160; m = n*sqrt(zeta^2/(1-zeta^2));
axis ([ -4 1 -4 4]);
plot (-m,n,'--'); % zeta
plot (-x,y,'rd');
saveas(fig,'Q2_Lead+Lag_rlocus.png');

fig = figure(6);
step(Gcl_lead,Gpf*Gcl_lead,Gcl_lead_lag,Gpf*Gcl_lead_lag)
legend('Lead','Pre-filter + Lead','Lead + Lag','Pre-filter + Lead + Lag'...
    ,'Location','southeast');
saveas(fig,'Q2_All_step_1.png');
stepinfo(Gpf*Gcl_lead_lag)

fig = figure(7);
margin(Gcl_lead_lag);
saveas(fig,'Q2_Lead+lag_margin.png');
BW_lead_lag = bandwidth(Gcl_lead_lag); % bandwidth

fig = figure(8);
margin(Gcl_lead_lag);
saveas(fig,'Q2_Lead+lag+prefilter_margin.png');
BW_lead_lag_prefilter = bandwidth(Gpf*Gcl_lead_lag); % bandwidth

%% Verifing Ka = 40
Ka = (z_lag/p_lag) * (z_lead/p_lead) * ( K/(a*b) );
ess_parabolic = a/Ka; 

%% New phase-lead compensator
z_lead2 = x; % given zero lead right below the desired pole
f = 180 + atand(y/(x-z_lead2))+(180-atand(y/(x-z_lag)))...
    +(180-atand(y/(x-z_lead))) - (180-atand(y/(x-p_lag)))-atand(y/(p_lead-x))...
    -2*(180-atand(y/x))-atand(y/(a-x))-atand(y/(b-x));
% p_lead2 = y/tand(180+f) - x;
p_lead2 = y/tand(f) + x;
Gc_lead2 = (s+z_lead2)/(s+p_lead2);
% Calculating the new gain K with the lead 2 compensator
K_new = sqrt((-x+p_lead2)^2+y^2)*sqrt((-x+p_lag)^2+y^2)*...
    sqrt((-x+p_lead)^2+y^2)*((-x)^2+y^2)*sqrt((-x+a)^2+y^2)*sqrt((-x+b)^2+y^2)/...
    ( sqrt((-x+z_lead)^2+y^2)*sqrt((-x+z_lag)^2+y^2)*sqrt((-x+z_lead)^2+y^2)*K );
Gol_lead2_lead_lag = K_new*Gc_lead2*Gol_lead_lag;
Gcl_lead2_lead_lag = feedback(Gol_lead2_lead_lag,1);
%% Verifing Ka = 40 
Ka_new = K_new*(z_lead2/p_lead2)*(z_lag/p_lag)*(z_lead/p_lead)*(K/(a*b));
ess_parabolic_new = a/Ka_new; 
%% Comparing the new compensator of the previous designs.
fig = figure(9);
rlocus(Gcl_lead2_lead_lag)
hold on;
% ploting the s1 and zeta in the rlocus
n = 0:1:160; m = n*sqrt(zeta^2/(1-zeta^2));
axis ([ -4 1 -4 4]);
plot (-m,n,'--'); % zeta
plot (-x,y,'rd');
saveas(fig,'Q2_Lead2+Lead+Lag_rlocus.png');

fig = figure(10);
% step(Gcl_lead2_lead_lag,Gpf*Gcl_lead2_lead_lag)
step(Gcl_lead,Gpf*Gcl_lead,Gcl_lead_lag,Gpf*Gcl_lead_lag,...
    Gcl_lead2_lead_lag,Gpf*Gcl_lead2_lead_lag)
legend('Lead','Pre-filter+Lead','Lead+Lag','Pre-filter+Lead+Lag'...
    ,'Lead 2+Lead+Lag','Pre-filter+Lead2+Lead+Lag','Location','southeast');
saveas(fig,'Q2_All_step_2.png');
stepinfo(Gpf*Gcl_lead2_lead_lag)

fig = figure(11);
margin(Gcl_lead2_lead_lag);
saveas(fig,'Q2_Lead2+lead+lag_margin.png');
BW_lead2_lead_lag = bandwidth(Gcl_lead2_lead_lag); % bandwidth

fig = figure(12);
margin(Gpf*Gcl_lead2_lead_lag);
saveas(fig,'Q2_Lead2+lead+lag+prefilter_margin.png');
BW_lead2_lead_lag_prefilter = bandwidth(Gpf*Gcl_lead2_lead_lag); % bandwidth

