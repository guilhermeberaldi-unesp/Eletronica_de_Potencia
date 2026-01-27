%% Retificador meia-onda + R // C (filtro capacitivo) + corrente no diodo i_D(t)
clear; clc; close all;

% ---- Parâmetros ----
f   = 60;
Vm  = 170;
R   = 5;
C   = 1470e-6;
Vd  = 0.7;
ncy = 3;
fs  = 200e3;

T  = 1/f;
dt = 1/fs;
t  = 0:dt:(ncy*T);
N  = numel(t);

v_s = Vm*sin(2*pi*f*t);

% ---- Simulação + detecção de condução ----
v_0    = zeros(1,N);
cond   = false(1,N);  % 1 quando o diodo conduz
for k = 2:N
    if v_s(k) > v_0(k-1) + Vd
        v_0(k)  = v_s(k) - Vd;    % carga
        cond(k) = true;
    else
        v_0(k)  = v_0(k-1)*exp(-dt/(R*C)); % descarga
    end
end

i_0 = v_0 / R;          % corrente na carga
v_D = v_s - v_0;        % tensão no diodo (anodo->catodo)

% ---- Corrente no capacitor e no diodo ----
dv0_dt = [0 diff(v_0)/dt];     % derivada numérica
i_C    = C * dv0_dt;           % i_C = C dv/dt

i_D = zeros(1,N);
i_D(cond) = i_0(cond) + i_C(cond);   % só existe quando conduz (modelo simples)
i_D(i_D<0) = 0;                      % segurança numérica (evita negativos por diff)

% ---- Paleta ----
statale.maincolor = [0   51  102]/255;
statale.lilla     = [120 0   80 ]/255;
statale.darkgreen = [0   70  40 ]/255;
statale.red       = [190 60  55 ]/255;
statale.yellow    = [200 155 20 ]/255;

plotColors = [
    statale.maincolor
    statale.lilla
    statale.darkgreen
    statale.red
    statale.yellow
];
set(groot,'defaultAxesColorOrder',plotColors);

LW = 3.0;
theta = 2*pi*f*t;

%% ---- i_D(t) FORÇADO (ilustrativo): pulsos todos iguais com pico ajustável ----
% Requer: t, f, ncy já definidos. Não mexe em v_s, v_0, i_0, v_D.
% Ajustes:
Ipk   = 100;        % <-- pico desejado (A)
Delta = 0.02*(1/f); % <-- meia-largura do pulso (s). Ex.: 0.01*T ... 0.05*T

T     = 1/f;
omega = 2*pi*f;

% instantes dos picos positivos (centro dos pulsos)
t_peaks = (pi/2 + 2*pi*(0:ceil(ncy)))/omega;
t_peaks = t_peaks(t_peaks >= t(1) & t_peaks <= t(end));

i_D = zeros(size(t));

for kk = 1:numel(t_peaks)
    tau = t - t_peaks(kk);
    m   = abs(tau) <= Delta;

    % pulso raised-cosine normalizado para ter máximo = 1 no centro
    g = zeros(size(t));
    g(m) = 0.5*(1 + cos(pi*tau(m)/Delta));  % g(0)=1, g(±Delta)=0

    i_D = i_D + Ipk*g;  % todos os pulsos com o mesmo pico Ipk
end



figure('Color','w');
plot(theta, v_s, 'Color', plotColors(1,:), 'LineWidth', LW); hold on;
plot(theta, v_0, 'Color', plotColors(2,:), 'LineWidth', LW);
plot(theta, i_0, 'Color', plotColors(3,:), 'LineWidth', LW);
plot(theta, v_D, 'Color', plotColors(4,:), 'LineWidth', LW);
plot(theta, i_D, 'Color', plotColors(5,:), 'LineWidth', LW);   % i_D(t) (yellow)
grid on;

xlabel('\theta (rad)');
ylabel('Tensão (V) / Corrente (A)');
title('Retificador monofásico de meia onda com carga resistiva e filtro capacitivo');

xlim([0, ncy*2*pi]);
xticks(0:pi:(ncy*2*pi));
xticklabels(compose('%g\\pi', (0:1:(ncy*2))));
xticklabels(strrep(xticklabels,'0\pi','0'));

lgd = legend('v_s(t)','v_0(t)','i_0(t)','v_D(t)','i_D(t)','Location','best');
lgd.FontSize   = 12;
lgd.FontWeight = 'bold';
lgd.Color      = 'none';

ax = gca;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontWeight = 'bold';

%% =========================
%  5 MP4s a partir do template .fig (sem precisar do script original)
%  Template: "meia onda rc.fig"
%  Resolução de render: 300 dpi (exportgraphics)
% =========================

templateFig = "meia onda rc.fig";
if ~isfile(templateFig)
    error("Não achei '%s' em %s", templateFig, pwd);
end

% Config vídeo
fps     = 30;
durEach = 10;                  % segundos por vídeo
speed = 2;  % 2x mais rápido
nFrames = max(60, round(durEach*fps/speed));
resDPI  = 300;                 % <-- pedido: 300 dpi
quality = 100;

names = ["v_s(t)","v_0(t)","i_0(t)","v_D(t)","i_D(t)"];
outFiles = ["rc_v1_vs.mp4", ...
            "rc_v2_v0.mp4", ...
            "rc_v3_i0.mp4", ...
            "rc_v4_vD.mp4", ...
            "rc_v5_iD.mp4"];

% ---------- 1) Lê os dados direto do .fig (uma vez) ----------
fig0 = openfig(templateFig, "new", "invisible");
drawnow;

% acha axes principal (evita axes da legenda)
axesList = findall(fig0, "Type","axes");
ax0 = [];
for a = axesList(:).'
    if isgraphics(a.XLabel) && isgraphics(a.YLabel)
        xl = string(a.XLabel.String);
        yl = string(a.YLabel.String);
        if contains(xl, "\theta") && contains(yl, "Tensão")
            ax0 = a; break;
        end
    end
end
if isempty(ax0), ax0 = axesList(1); end

plines0 = findall(ax0, "Type","line");
dn0 = string(get(plines0,"DisplayName"));

hs0 = gobjects(1,5);
for j = 1:5
    idx = find(dn0 == names(j), 1, "first");
    if isempty(idx)
        close(fig0);
        error("Não achei a curva '%s' no template (DisplayName).", names(j));
    end
    hs0(j) = plines0(idx);
end

% dados vindos do .fig
X1 = hs0(1).XData(:);
Y  = zeros(numel(X1), 5);
for j = 1:5
    Y(:,j) = hs0(j).YData(:);
end
N = numel(X1);

% guarda tamanho original da figura (pra ficar igual)
set(fig0,'Units','pixels');
pos0 = get(fig0,'Position');      % [x y w h]
pos0(3) = pos0(3) + mod(pos0(3),2);  % w par (H.264)
pos0(4) = pos0(4) + mod(pos0(4),2);  % h par
close(fig0);

reveal = @(y,idx) [y(1:idx); nan(N-idx,1)];

% ---------- 2) Gera 5 vídeos ----------
for vid = 1:5
    fig = openfig(templateFig, "new", "visible");
    drawnow;

    % trava tamanho/layout conforme template (ajustado para w/h pares)
    set(fig,'Units','pixels');
    set(fig,'Resize','off');
    set(fig,'Position',pos0);
    drawnow;

    % acha axes principal
    axesList = findall(fig, "Type","axes");
    ax = [];
    for a = axesList(:).'
        if isgraphics(a.XLabel) && isgraphics(a.YLabel)
            xl = string(a.XLabel.String);
            yl = string(a.YLabel.String);
            if contains(xl, "\theta") && contains(yl, "Tensão")
                ax = a; break;
            end
        end
    end
    if isempty(ax), ax = axesList(1); end

    % acha linhas por DisplayName
    plines = findall(ax, "Type","line");
    dn = string(get(plines,"DisplayName"));
    

    hs = gobjects(1,5);
    for j = 1:5
        idx = find(dn == names(j), 1, "first");
        if isempty(idx)
            close(fig);
            error("Não achei a curva '%s' no template (DisplayName).", names(j));
        end
        hs(j) = plines(idx);
    end

    % fixa XData (preserva estilo)
    set(hs, "XData", X1);

    % estado inicial do vídeo:
    for j = 1:5
        if j < vid
            set(hs(j), "YData", Y(:,j));       % anteriores completas
        else
            set(hs(j), "YData", nan(N,1));     % atual (animará) + futuras vazias
        end
    end

    % ---- FIXA ESCALA (Y e X) para não variar durante a animação ----
yMin = min(Y(:));
yMax = max(Y(:));
pad  = 0.05*(yMax - yMin + eps);
set(ax, 'YLim', [yMin-pad, yMax+pad], 'YLimMode', 'manual');

% (recomendado também) fixa XLim
set(ax, 'XLim', [X1(1) X1(end)], 'XLimMode', 'manual');

    drawnow;

    % VideoWriter
    vw = VideoWriter(outFiles(vid), "MPEG-4");
    vw.FrameRate = fps;
    vw.Quality   = quality;
    open(vw);
    cleanupObj = onCleanup(@() close(vw));  % evita MP4 0 KB se der erro

    % temp png único por vídeo
    tmpPng = fullfile(tempdir, "frame_tmp_" + string(feature("getpid")) + "_v" + vid + ".png");

    % frame teste (define HxW de referência)
    exportgraphics(fig, tmpPng, 'Resolution', resDPI);
    im0 = imread(tmpPng);
    if size(im0,3) == 1, im0 = repmat(im0,1,1,3); end
    H = size(im0,1); W = size(im0,2);

    for k = 1:nFrames
        frac = (k-1)/(nFrames-1);
        idx  = max(2, round(frac*(N-1)) + 1);

        set(hs(vid), "YData", reveal(Y(:,vid), idx));
        drawnow;

        exportgraphics(fig, tmpPng, 'Resolution', resDPI);
        im = imread(tmpPng);
        if size(im,3) == 1, im = repmat(im,1,1,3); end

        % normaliza tamanho fixo (crop/pad sem toolbox)
        im = im(1:min(end,H), 1:min(end,W), :);
        dh = H - size(im,1);
        dw = W - size(im,2);
        if dh > 0 || dw > 0
            im2 = uint8(255*ones(H, W, 3)); % fundo branco
            im2(1:size(im,1), 1:size(im,2), :) = im;
            im = im2;
        end
        im = im(1:end-mod(size(im,1),2), 1:end-mod(size(im,2),2), :);

        writeVideo(vw, im);
    end

    if exist(tmpPng,'file'), delete(tmpPng); end
    clear cleanupObj;  % fecha vw
    close(fig);

    disp("Salvou: " + fullfile(pwd, outFiles(vid)));
end
