%% Set up axes
w = logspace(-2,2);
t = linspace(0,40);

% Define the plant
nG = 9.8*[1 -.5 6.3];
dG = conv([1 .6565],[1 -.2366 .1493]);

G = tf(nG,dG);
bode(G,w);



%% Crossover Frequency

wc = 1; %desired crossover frequency
[mag, ph] = bode(G,wc);
K0 = 1/mag;
bode(G,K0*G,w);


%% Lead

a = sqrt(10);

% Make a first order lead
K1 = tf(a*[1 wc/a],[1 a*wc]);

% Check, should give mag = 1, ph approx 55 degrees
[mag,ph] = bode(K1,wc);
bode(K1,w);

% now form the full compensator
K = K0*K1*K1;
bode(K,w);

bode(G*K0,G*K,w);

%% Return ratio
L = K*G;

%check, should give mag=1, ph approx -150 degrees
[mag,ph] = bode(L,wc);


%% Complementary sensitivity function "T"
T = feedback(L,1);

% closed loop poles
pole(T);

% or in a nicer form
% damp(T);

% look at the step respons
step(T,t);

% and its frequency response
bode(T,w);

%% Closed loop transfer function g/1 +gk
CL = feedback(G,K);

% look at the step response
step(CL,t);

dcgain = bode(CL,0);

% form the precompensator to make closed-loop DC gain = 1
Ka = 1/dcgain;

% step response
step(Ka*CL,t);

% and its frequency response
bode(Ka*CL,w);

%preturbaton
nG1 = 9*[1 -0.5 6.3];

G1 = tf(nG1,dG);
CL1 = feedback(G1,K);

step(Ka*CL1,t);


%% Steady state error
%try adding some phaes lag
K2 = tf([1 .1 ],[1 .002]);
bode(K2,w); %looks ok

K = K2*K0*K1*K1;
CL = feedback(G,K);

dcgain = bode(CL,0); 
Ka = 1/dcgain;
step(Ka*CL); %not so good. look at a bode plot to see why
bode(Ka*CL);

%% Try putting the lag compensator in the forward path instead

Kb = K0*K1*K1;
Kc = K2;

CL = feedback(G*Kc,Kb);
dcgain = bode(CL,0); Ka = 1/dcgain;
Ka = 0.0021;
CL = feedback(G*Kc,Kb)*Ka;
bode(CL);

%% How bout the step response?
step(CL,t);

% better but slow
% Modify Ka directly to speed up response
% Or even better, design the desired tracking response directly

% CLref must include unstable zeros, for internal stability.
tzero(CL);
roots(nG);
CLref = tf(nG/nG(3), poly([-1 -1 -1 -1]));
step(CLref,t);

Ka = CLref/CL;
tzero(Ka);
pole(Ka);

Ka = minreal(Ka,1e-3); %find minimal realization
pole(Ka);
CL = feedback(G*Kc,Kb)*Ka/dcgain;
step(CL,t);

% that's better
% now check robustness

CL1 = feedback(G1*Kc,Kb)*Ka/dcgain;
step(CL1,t);
%% *****Increasing line weight****
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
grid on
%************************************