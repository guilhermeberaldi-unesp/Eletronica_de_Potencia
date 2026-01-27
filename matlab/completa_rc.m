%% Retificador monofásico de ONDA COMPLETA (ponte) + carga R + filtro capacitivo
% Plota: v_s(t), v_rec(t), v_0(t), i_0(t) e i_D(t) (FORÇADO matematicamente)
% Cores (ordem): maincolor, lilla, darkgreen, red, yellow
% Eixo x em theta (rad) no recorte, começando em 0 com ticks em múltiplos de pi

clear; clc; close all;

% ---- Parâmetros ----
f   = 60;           % Hz
Vm  = 170;          % Vp
R   = 7;           % ohms
C   = 670e-6;       % F
Vd  = 0.7;          % V (queda por diodo) -> ponte: 2 diodos em condução
ncy = 12;           % ciclos simulados
fs  = 200e3;        % Hz

% ---- Tempo ----
T  = 1/f;
dt = 1/fs;
t  = 0:dt:(ncy*T);
N  = numel(t);

% ---- Fonte ----
omega = 2*pi*f;
theta = omega*t;
v_s   = Vm*sin(theta);

% ---- Retificado da ponte (disponível na saída antes do capacitor) ----
v_rec = abs(v_s) - 2*Vd;
v_rec(v_rec < 0) = 0;

% ---- Simulação do filtro capacitivo (modelo clamp/descarga) ----
% Se v_rec > v0_prev => carrega: v0 = v_rec
% Senão => descarga em R: v0(k)=v0(k-1)*exp(-dt/(RC))
v_0 = zeros(1,N);
for k = 2:N
    if v_rec(k) > v_0(k-1)
        v_0(k) = v_rec(k);
    else
        v_0(k) = v_0(k-1)*exp(-dt/(R*C));
    end
end

% ---- Corrente na carga ----
i_0 = v_0 / R;

%% ---- Corrente no diodo FORÇADA (matemática), igual ao meia-onda ----
% Para onda completa, há recarga a cada meio-período (T/2),
% então usamos picos de |sin| em theta = pi/2 + k*pi.
Ipk   = 150;            % <-- pico desejado do pulso (A)
Delta = 0.02*(T/2);    % <-- meia-largura do pulso (s)  (0.01*(T/2) ... 0.05*(T/2))

% picos de |sin| (ponte completa): theta = pi/2 + k*pi
t_peaks = (pi/2 + pi*(0:ceil(2*ncy)))/omega;
t_peaks = t_peaks(t_peaks >= t(1) & t_peaks <= t(end));

i_D = zeros(size(t));

for kk = 1:numel(t_peaks)
    t0  = t_peaks(kk);
    tau = t - t0;

    m = abs(tau) <= Delta;

    % pulso raised-cosine (máximo 1 no centro)
    g = zeros(size(t));
    g(m) = 0.5*(1 + cos(pi*tau(m)/Delta));   % g(0)=1

    i_D = i_D + Ipk*g;   % todos os pulsos com o mesmo pico
end

% ---- Paleta (RGB 0..1) ----
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

%% ---- ESCOLHA DO INTERVALO DE CICLOS PARA PLOTAR ----
plotCycles = [6 7];     % [ciclo_inicial ciclo_final]
c1 = plotCycles(1);
c2 = plotCycles(2);

t1 = (c1-1)*T;
t2 = c2*T;

idx = (t >= t1) & (t <= t2);

thetaP = theta(idx);
thetaRel = thetaP - thetaP(1);

vsP   = v_s(idx);
vrecP = v_rec(idx);
v0P   = v_0(idx);
i0P   = i_0(idx);
iDP   = i_D(idx);

%% ---- Plot ----
figure('Color','w');
plot(thetaRel, vsP,   'Color', plotColors(1,:), 'LineWidth', LW); hold on; % maincolor
plot(thetaRel, v0P,   'Color', plotColors(2,:), 'LineWidth', LW);          % lilla
plot(thetaRel, i0P,   'Color', plotColors(3,:), 'LineWidth', LW);          % darkgreen
plot(thetaRel, vrecP, 'Color', plotColors(4,:), 'LineWidth', LW);          % red
plot(thetaRel, iDP,   'Color', plotColors(5,:), 'LineWidth', LW);          % yellow
grid on;

xlabel('\theta (rad)');
ylabel('Tensão (V) / Corrente (A)');
title('Retificador monofásico de onda completa (ponte) com carga R e filtro capacitivo');

% ---- Ticks em múltiplos de pi a partir de 0 ----
xlim([0, thetaRel(end)]);
k2 = floor(thetaRel(end)/pi);
xt = (0:k2)*pi;
xticks(xt);

m = 0:k2;
labs = strings(size(m));
for kk = 1:numel(m)
    if m(kk) == 0
        labs(kk) = "0";
    elseif m(kk) == 1
        labs(kk) = "\pi";
    else
        labs(kk) = sprintf('%g\\pi', m(kk));
    end
end
xticklabels(labs);

% ---- Legenda + fontes ----
lgd = legend('v_s(t)','v_0(t)','i_0(t)','v_{rec}(t)','i_D(t) (forçado)', 'Location','best');
lgd.FontSize   = 12;
lgd.FontWeight = 'bold';
lgd.Color      = 'none';

ax = gca;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontWeight = 'bold';


%% =========================
%  5 MP4s a partir do template .fig (sem precisar do script original)
%  Template: "completa rc.fig"
%  Resolução de render: 300 dpi (exportgraphics)
%  2x mais rápido (metade dos frames)
%  V1: anima v_s
%  V2: v_s fixa + anima v_0
%  V3: v_s,v_0 fixas + anima i_0
%  V4: v_s,v_0,i_0 fixas + anima v_{ret}(t)
%  V5: v_s,v_0,i_0,v_{ret} fixas + anima i_D
%
%  OBS: para o HxW ficar coerente com a FIG:
%  - usamos OuterPosition do template
%  - fixamos OuterPosition antes do export
%  - normalizamos cada frame para HxW (crop/pad)
% =========================

templateFig = "completa rc.fig";
if ~isfile(templateFig)
    error("Não achei '%s' em %s", templateFig, pwd);
end

% Config vídeo
fps     = 30;
durEach = 10;        % segundos por vídeo (antes do "speed")
speed   = 2;         % 2x mais rápido
nFrames = max(60, round(durEach*fps/speed));
resDPI  = 300;       % <-- 300 dpi
quality = 100;

names = ["v_s(t)","v_0(t)","i_0(t)","v_{ret}(t)","i_D(t)"];
outFiles = ["comp_rc_v1_vs.mp4", ...
            "comp_rc_v2_v0.mp4", ...
            "comp_rc_v3_i0.mp4", ...
            "comp_rc_v4_vret.mp4", ...
            "comp_rc_v5_iD.mp4"];

% ---------- 1) Lê os dados e o OuterPosition direto do .fig (uma vez) ----------
fig0 = openfig(templateFig, "new", "invisible");
drawnow;

set(fig0,'Units','pixels');
outer0 = get(fig0,'OuterPosition');          % [x y w h]
outer0(3) = outer0(3) + mod(outer0(3),2);    % w par (H.264)
outer0(4) = outer0(4) + mod(outer0(4),2);    % h par (H.264)

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

% fixa limites globais (evita autoscale enquanto anima)
yMin = min(Y(:)); yMax = max(Y(:));
pad  = 0.05*(yMax - yMin + eps);
YLimFixed = [yMin-pad, yMax+pad];

close(fig0);

reveal = @(y,idx) [y(1:idx); nan(N-idx,1)];

% ---------- 2) Gera 5 vídeos ----------
for vid = 1:5
    fig = openfig(templateFig, "new", "visible");
    drawnow;

    % trava tamanho/layout exatamente pelo OuterPosition do template
    set(fig,'Units','pixels');
    set(fig,'Resize','off');
    set(fig,'OuterPosition', outer0);
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

    % fixa escala desde o começo
    set(ax, 'XLim', [X1(1) X1(end)], 'XLimMode','manual');
    set(ax, 'YLim', YLimFixed,       'YLimMode','manual');

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
    drawnow;

    % VideoWriter
    vw = VideoWriter(outFiles(vid), "MPEG-4");
    vw.FrameRate = fps;
    vw.Quality   = quality;
    open(vw);
    cleanupObj = onCleanup(@() close(vw));  % evita MP4 0 KB se der erro

    % temp png único por vídeo
    tmpPng = fullfile(tempdir, "frame_tmp_" + string(feature("getpid")) + "_v" + vid + ".png");

    % frame teste (define HxW de referência a partir do export)
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
