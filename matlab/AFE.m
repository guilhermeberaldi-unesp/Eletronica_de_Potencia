%% AFE Trifásico (Retificador Ativo) - Modelo médio com controle dq
% Plots: v0 (Vdc), i0 (Idc), vs_abc e i_abc de entrada
% Modelo: rede trifásica -> indutor Lf/Rf -> conversor (tensão controlada v_conv_abc) -> barramento CC (Cdc + carga)
%
% Autor: ChatGPT
% -------------------------------------------------------------

clear; clc; close all;

%% ====== Parâmetros elétricos ======
f    = 60;                 % Hz
w    = 2*pi*f;
Vllr = 380;                % V RMS linha-linha (rede)
Vphr = Vllr/sqrt(3);       % V RMS fase-neutro
Vphm = sqrt(2)*Vphr;       % V pico fase-neutro

Lf   = 2.5e-3;             % H (indutor por fase)
Rf   = 0.15;               % ohm (resistência série por fase)

Cdc  = 10e-3;               % F (capacitor do barramento CC)
Rload= 50;                % ohm (carga equivalente no barramento)  -> i0 = Vdc/Rload

Vdc_ref = 700;             % V referência do barramento (>= ~1.35*Vllr para folga)

%% ====== Controle ======
% PLL ideal: theta = w*t (rede perfeita)
% Malha de corrente (dq): PI
Kpi = 8;                   % ganho proporcional corrente
Kii = 2000;                % ganho integral corrente

% Malha de tensão do barramento: PI gera i_d* (referência de corrente ativa)
Kpv = 0.2;                 % ganho proporcional tensão
Kiv = 2;                  % ganho integral tensão

% Fator de potência ~1: i_q* = 0
iq_ref = 0;

% Limite de corrente (proteção)
Imax = 100;                 % A (pico em dq)

%% ====== Simulação ======
Ts   = 5e-6;               % passo (50 us)
Tsim = 5;               % s
t    = 0:Ts:Tsim;
N    = numel(t);

% Estados
iabc  = zeros(3,N);        % [ia; ib; ic]
Vdc   = zeros(1,N);        % v0(t)
Vdc(1)= 400;               % condição inicial

% Integradores PI
xi_d = 0; xi_q = 0;        % integradores das correntes
xv   = 0;                  % integrador da tensão do barramento

% Logs extras
vsabc   = zeros(3,N);
idq_log = zeros(2,N);
vdqcmd  = zeros(2,N);
idref_log = zeros(1,N);

%% ====== Funções auxiliares (Clarke/Park) ======
clarke = @(abc) (2/3)*[ 1      -1/2     -1/2;
                       0  sqrt(3)/2 -sqrt(3)/2 ]*abc;

invclarke = @(ab) [ 1,            0;
                   -1/2,  sqrt(3)/2;
                   -1/2, -sqrt(3)/2 ]*ab;

park = @(ab,th) [ cos(th)  sin(th);
                 -sin(th)  cos(th) ]*ab;

invpark = @(dq,th) [ cos(th) -sin(th);
                     sin(th)  cos(th) ]*dq;

%% ====== Loop de simulação ======
for k = 1:N-1

    th = w*t(k);

    % Rede trifásica (fase-neutro)
    vsa = Vphm*cos(th);
    vsb = Vphm*cos(th - 2*pi/3);
    vsc = Vphm*cos(th + 2*pi/3);
    vs  = [vsa; vsb; vsc];
    vsabc(:,k) = vs;

    % Corrente de carga no barramento (i0)
    i0 = Vdc(k)/Rload;

    % (aprox.) potência no barramento = Vdc*i0
    % i_d* via PI de tensão (controle de energia)
    ev = Vdc_ref - Vdc(k);
    xv = xv + Ts*ev;
    id_ref = Kpv*ev + Kiv*xv;

    % Limitar referência
    id_ref = max(min(id_ref, Imax), -Imax);
    idref_log(k) = id_ref;

    % Medições: iabc -> i_alpha_beta -> i_dq
    i_ab = clarke(iabc(:,k));
    i_dq = park(i_ab, th);
    idq_log(:,k) = i_dq;

    % Erros de corrente
    ed = id_ref - i_dq(1);
    eq = iq_ref - i_dq(2);

    % PI corrente
    xi_d = xi_d + Ts*ed;
    xi_q = xi_q + Ts*eq;

    % Termo de desacoplamento + queda no indutor (feedforward)
    % Dinâmica: L di/dt = v_s - v_conv - R i  (em abc)
    % Em dq (com rede ideal): L did/dt = v_sd - v_cd - R id + wL iq
    %                          L diq/dt = v_sq - v_cq - R iq - wL id
    % Logo, comando v_cd = v_sd - R id + wL iq - u_d
    %       comando v_cq = v_sq - R iq - wL id - u_q
    % onde u_d = PI(ed), u_q = PI(eq)
    v_ab = clarke(vs);
    v_dq = park(v_ab, th);
    v_sd = v_dq(1);
    v_sq = v_dq(2);

    u_d = Kpi*ed + Kii*xi_d;
    u_q = Kpi*eq + Kii*xi_q;

    v_cd = v_sd - Rf*i_dq(1) + w*Lf*i_dq(2) - u_d;
    v_cq = v_sq - Rf*i_dq(2) - w*Lf*i_dq(1) - u_q;

    vdqcmd(:,k) = [v_cd; v_cq];

    % Volta para abc: v_conv_abc
    v_cab = invpark([v_cd; v_cq], th);
    v_cabc = invclarke(v_cab);

    % Dinâmica das correntes em abc: di/dt = (vs - vconv - R i)/L
    di = (vs - v_cabc - Rf*iabc(:,k))/Lf;
    iabc(:,k+1) = iabc(:,k) + Ts*di;

    % Potência AC (instantânea) e dinâmica do barramento
    Pac(k) = vs.' * iabc(:,k);                 % (aprox) potência trifásica instantânea
    dVdc = (Pac(k) - Vdc(k)*i0) / (Cdc*max(Vdc(k),1));  
    % Observação: relação energia: d/dt(1/2 C V^2)= Pac - Pload -> C V dV/dt = Pac - Pload
    Pload(k) = (Vdc(k)^2)/Rload;
    Vdc(k+1) = Vdc(k) + Ts*dVdc;
    Vdc_safe = max(Vdc(k), 50);      % evita divisão por ~0
    dVdc = (Pac(k) - Pload(k)) / (Cdc * Vdc_safe);

    Vdc(k+1) = Vdc(k) + Ts*dVdc;

end

% guardar último ponto
vsabc(:,N) = vsabc(:,N-1);
idq_log(:,N) = idq_log(:,N-1);

%% ====== Sinais de saída ======
v0 = Vdc;              % tensão CC
i0 = Vdc./Rload;       % corrente CC
ia = iabc(1,:); ib = iabc(2,:); ic = iabc(3,:);

%% ====== Plots (mesmo estilo Statale que estávamos usando) ======

% ---- Paleta (RGB 0..1) ----
statale.maincolor = [0   51  102]/255;
statale.lilla     = [120 0   80 ]/255;
statale.darkgreen = [0   70  40 ]/255;
statale.red       = [190 60  55 ]/255;
statale.yellow    = [200 155 20 ]/255;
statale.grey      = [50  50  50 ]/255;

plotColors = [
    statale.maincolor
    statale.lilla
    statale.darkgreen
    statale.red
    statale.yellow
    statale.grey
];
set(groot,'defaultAxesColorOrder',plotColors);

LW3 = 3.0;   % principais
LW2 = 2.0;   % secundárias

FSleg = 12;

figure('Color','w','Position',[80 80 1100 750]);
tiledlayout(2,1,"Padding","compact","TileSpacing","compact");

% (1) Barramento CC
nexttile;
plot(t, v0, 'LineWidth', LW3, 'Color', plotColors(1,:)); grid on;
ylabel('v_0 = V_{dc} (V)','FontWeight','bold');
title('AFE trifásico - Barramento CC');

yyaxis right;
plot(t, i0, 'LineWidth', LW3, 'Color', plotColors(6,:)); grid on;
ylabel('i_0 = I_{dc} (A)','FontWeight','bold');

lgd = legend('v_0','i_0','Location','best');
lgd.FontSize = FSleg; lgd.FontWeight = 'bold'; lgd.Color = 'none';

ax = gca;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontWeight = 'bold';

% % (2) Tensões de entrada
% nexttile;
% plot(t, vsabc(1,:), 'LineWidth', LW3, 'Color', plotColors(1,:)); hold on;
% plot(t, vsabc(2,:), 'LineWidth', LW3, 'Color', plotColors(2,:));
% plot(t, vsabc(3,:), 'LineWidth', LW3, 'Color', plotColors(3,:));
% grid on; 
% ylabel('v_s (V)','FontWeight','bold');
% title('Tensões de entrada (v_{sa}, v_{sb}, v_{sc})');
% 
% lgd = legend('v_{sa}','v_{sb}','v_{sc}','Location','best');
% lgd.FontSize = FSleg; lgd.FontWeight = 'bold'; lgd.Color = 'none';
% 
% ax = gca;
% ax.XLabel.FontWeight = 'bold';
% ax.YLabel.FontWeight = 'bold';

% (3) Correntes de entrada
nexttile;
plot(t, ia, 'LineWidth', LW3, 'Color', plotColors(3,:)); hold on;
plot(t, ib, 'LineWidth', LW3, 'Color', plotColors(4,:));
plot(t, ic, 'LineWidth', LW3, 'Color', plotColors(5,:));
grid on; 
xlabel('Tempo (s)','FontWeight','bold'); 
ylabel('i (A)','FontWeight','bold');
title('Correntes de entrada (i_a, i_b, i_c)');

lgd = legend('i_a','i_b','i_c','Location','best');
lgd.FontSize = FSleg; lgd.FontWeight = 'bold'; lgd.Color = 'none';

ax = gca;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontWeight = 'bold';

%% ====== (Opcional) dq e referência (mesmo estilo) ======
figure('Color','w','Position',[150 150 1100 420]);
plot(t, idq_log(1,:), 'LineWidth', LW3, 'Color', plotColors(1,:)); hold on;
plot(t, idq_log(2,:), 'LineWidth', LW3, 'Color', plotColors(2,:));
plot(t, idref_log, '--', 'LineWidth', LW2, 'Color', plotColors(5,:));
grid on; 
xlabel('Tempo (s)','FontWeight','bold'); 
ylabel('Corrente (A)','FontWeight','bold');
title('Correntes em dq (controle)');

lgd = legend('i_d','i_q','i_d^*','Location','best');
lgd.FontSize = FSleg; lgd.FontWeight = 'bold'; lgd.Color = 'none';

ax = gca;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontWeight = 'bold';
