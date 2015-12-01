clear;
% for alpha = 0.045:0.010:0.055
for beta = 0:0
    
    % Ka = 1;
    % Kb = tf([1.3 1],[1 1.3]);
    % Kc = 1;
    wc = 1;
%     Ka = 1;
    
%% Lead compensator
    lead = 55;
    lead_coeff = tand(0.5*(90+lead));
    Kb = lead_coeff*tf([1 wc/lead_coeff],[1 wc*lead_coeff]) ;
    
%% Lag compensator  
    wc = 5/10;
    lag = 10;
    lag_coeff = tand(0.5*(90+lag));
    Kc =  1.2*(1/lag_coeff)*tf([1 wc*lag_coeff],[1 wc/lag_coeff]);
%     Kc = 1;
% %     Kc = tf([1],[1  0]);
%     Kc = 1;
%     Kc = 7;
%     Kc = 1/(20.004*(1/1.0008 -1/sqrt(10))); %0.0732
%       Kc = (1/0.0732)*tf([1  wc*0.0732],[1  wc/0.0732]);
% Kc = 0.0732;

%% Integrator
% wc = 10;
% T = 1/wc;
% Kc = tf([2],[T^2 2*T 1] );

%% Plant
    alpha = 0.05;
    % beta = 0;
    
    % first form the plant as a sum of G1 and G2
    G1 = tf([1],[1 -0.1 alpha]);
    G2 = tf([beta 0.1],[1 3 25]);
    
    G = G1+G2;
    
    OL = G*Kc*Kb;
    CL = feedback(G*Kc,Kb); %CLTF from output of Ka to y
    
    %% Precompensator
    % To meet requirements of 2
    dcgain = bode(CL,0);
    Ka  =1/dcgain;
%     Ka = 1;
%     Ka = 1/1.0008;
    
    CLref = CL*Ka; %this is CLTF from r to y
    
    hold on
    figure(1); margin(OL);legend('\beta = -1, \alpha = 0.05','\beta = 0, \alpha = 0.05','\beta = 1, \alpha = 0.05'); %gain and ph margins
    figure(2); step(CLref); legend('\beta = -1, \alpha = 0.05','\beta = 0, \alpha = 0.05','\beta = 1, \alpha = 0.05');%step response
    hold on
    figure(3);nyquist(OL);legend('\beta = -1, \alpha = 0.05','\beta = 0, \alpha = 0.05','\beta = 1, \alpha = 0.05');
    
    %% Internal Stability
    
    % Show that magnitude of T remains below 1/delta
    T = tf([Kb*G],[Kb*G+1]);
    delta = tf([conv([beta +0.1],[1 0.1 + 0.05])],[1 3 25 0]);
    
    %     figure(4);bode(T,1/delta);legend('T','1/\Delta s')
    
    % Check S etc
    
    
    %% Precompensator
    %     %To meet requirements of 2
    %
    %     dcgain = bode(CL,0);
    %     Ka  =1/dcgain;
    %     hold on
    %     figure(5);step(Ka*CL);legend('\beta = -1, \alpha = 0.05','\beta = 0, \alpha = 0.05','\beta = 1, \alpha = 0.05');
    %
    
    % Looking at final value (SSE)
    t = 0:0.1:50;
    u = t;
    [y,t,x] = lsim(CLref,u,t);
    figure(6)
    hold on;
    plot(t,y,'y',t,u,'m')
    xlabel('Time (sec)')
    ylabel('Amplitude')
    title('Input-purple, Output-yellow')
    
    % Steady state error.
    % try adding some phase lag
    
    %eg 20 degrees--> alpha = 1.428
    %     Kc = (1/1.428)*tf([1 1.428],[1 1/1.428]);
    
    
    figure(7);pzmap(OL);hold on
    
%     %% *****Increasing line weight****
%     h = findobj(gcf,'type','line');
%     set(h,'linewidth',2);
%     grid on
%     %************************************
end
% end