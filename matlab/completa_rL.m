%% Retificador monofásico de ponte completa + carga RL + diodo de roda livre (Df)
% Plota: v_s(t), v_0(t), i_0(t) e i_D(t) (corrente no "ramo retificador" da ponte)
% Observação: i_D(t) aqui representa a corrente que passa pelo retificador (ponte) quando a fonte alimenta a carga.
% Durante a roda livre (Df), i_D = 0.
% Cores (ordem): maincolor, lilla, darkgreen, red, yellow
% Eixo x em theta (rad) com ticks em múltiplos de pi (começando em 0 no recorte)

clear; clc; close all;

% ---- Parâmetros ----
f   = 60;          % Hz
Vm  = 170;         % Vp
R   = 2;           % ohms
L   = 5e-5;       % H
Vd  = 0.7;         % V (queda por diodo)
Vdf = 0.7;         % V (diodo de roda livre)
ncy = 30;          % ciclos simulados
fs  = 200e3;       % Hz

% ---- Tempo ----
T  = 1/f;
dt = 1/fs;
t  = 0:dt:(ncy*T);
N  = numel(t);

% ---- Fonte ----
omega = 2*pi*f;
theta = omega*t;
v_s   = Vm*sin(theta);

% ---- Inicialização ----
i_0 = zeros(1,N);   % corrente na carga (R-L)
v_0 = zeros(1,N);   % tensão na carga
i_D = zeros(1,N);   % corrente no retificador (ponte) quando conduz

% ---- Simulação por estados ----
% Ponte completa: quando alimenta pela rede, a tensão na carga é ~|v_s| - 2Vd
% Roda livre: v0 ~ -Vdf (aprox 0 para ideal), corrente decai em R-L
for k = 2:N
    i_prev = i_0(k-1);

    vrec = abs(v_s(k)) - 2*Vd;   % tensão disponível na saída da ponte
    if vrec < 0
        vrec = 0;
    end

    if (vrec > -Vdf) && (vrec > 0)  % ponte consegue alimentar
        v_0(k) = vrec;

        di   = (v_0(k) - R*i_prev)/L;
        iNew = i_prev + di*dt;

        if iNew < 0, iNew = 0; end
        i_0(k) = iNew;
        i_D(k) = iNew;            % corrente vem da ponte

        if i_0(k) == 0
            v_0(k) = 0;
            i_D(k) = 0;
        end

    elseif i_prev > 0
        % roda livre
        v_0(k) = -Vdf;

        di   = (v_0(k) - R*i_prev)/L;
        iNew = i_prev + di*dt;

        if iNew < 0, iNew = 0; end
        i_0(k) = iNew;
        i_D(k) = 0;

        if i_0(k) == 0
            v_0(k) = 0;
        end

    else
        v_0(k) = 0;
        i_0(k) = 0;
        i_D(k) = 0;
    end
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
plotCycles = [10 12];
c1 = plotCycles(1);
c2 = plotCycles(2);

t1 = (c1-1)*T;
t2 = c2*T;

idx = (t >= t1) & (t <= t2);

thetaP = theta(idx);
vsP    = v_s(idx);
v0P    = v_0(idx);
i0P    = i_0(idx);
iDP    = i_D(idx);

% eixo x relativo no recorte (começa em 0)
thetaRel = thetaP - thetaP(1);

%% ---- Plot ----
figure('Color','w');
plot(thetaRel, vsP, 'Color', plotColors(1,:), 'LineWidth', LW); hold on;
plot(thetaRel, v0P, 'Color', plotColors(2,:), 'LineWidth', LW);
plot(thetaRel, i0P, 'Color', plotColors(3,:), 'LineWidth', LW);
plot(thetaRel, iDP, 'Color', plotColors(4,:), 'LineWidth', LW);
grid on;

xlabel('\theta (rad)');
ylabel('Tensão (V) / Corrente (A)');
title('Retificador monofásico de ponte completa com carga RL e diodo de roda livre');

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


%% ---- Linha média (aprox.): V0_avg ≈ 2*Vmax/pi se Vmax >> Vpol ----
Vmax = Vm;        % pico da fonte
Vpol = 2*Vd;      % "tensão de polarização" da ponte (2 diodos). Ajuste se quiser outro critério.

Vavg_approx = 2*Vmax/pi;

hAvg = yline(Vavg_approx, '--', 'LineWidth', 2.5);
hAvg.Color = plotColors(5,:);  % yellow

txt = '$V_{0,\mathrm{avg}} \approx \frac{2V_{\max}}{\pi}\ \ \mathrm{se}\ \ V_{\max} \gg V_{\mathrm{pol}}$';
text(0.02*thetaRel(end), 0.92*max(v0P), txt, ...
    'Interpreter','latex', 'FontWeight','bold', 'FontSize', 12, ...
    'VerticalAlignment','top');

lgd = legend('v_s(t)','v_0(t)','i_0(t)','i_D(t)','V_{med}','Location','best');
lgd.FontSize   = 12;
lgd.FontWeight = 'bold';
lgd.Color      = 'none';

ax = gca;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontWeight = 'bold';

%% =========================
%  3 MP4s a partir do template .fig (sem precisar do script original)
%  Template: "onda completa rL.fig"
%  Resolução de render: 300 dpi (exportgraphics)
%  2x mais rápido (metade dos frames)
%  V1: anima v_s
%  V2: v_s fixa + anima v_0
%  V3: v_s,v_0 fixas + anima i_0
% =========================

templateFig = "onda completa rL.fig";
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

names = ["v_s(t)","v_0(t)","i_0(t)"];
outFiles = ["ponte_rL_v1_vs.mp4", ...
            "ponte_rL_v2_v0.mp4", ...
            "ponte_rL_v3_i0.mp4"];

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

hs0 = gobjects(1,3);
for j = 1:3
    idx = find(dn0 == names(j), 1, "first");
    if isempty(idx)
        close(fig0);
        error("Não achei a curva '%s' no template (DisplayName).", names(j));
    end
    hs0(j) = plines0(idx);
end

% dados vindos do .fig
X1 = hs0(1).XData(:);
Y  = zeros(numel(X1), 3);
for j = 1:3
    Y(:,j) = hs0(j).YData(:);
end
N = numel(X1);

% fixa limites globais (evita autoscale enquanto anima)
yMin = min(Y(:)); yMax = max(Y(:));
pad  = 0.05*(yMax - yMin + eps);
YLimFixed = [yMin-pad, yMax+pad];

% guarda tamanho original da figura (pra ficar igual)
set(fig0,'Units','pixels');
pos0 = get(fig0,'Position');          % [x y w h]
pos0(3) = pos0(3) + mod(pos0(3),2);   % w par (H.264)
pos0(4) = pos0(4) + mod(pos0(4),2);   % h par
close(fig0);

reveal = @(y,idx) [y(1:idx); nan(N-idx,1)];

% ---------- 2) Gera 3 vídeos ----------
for vid = 1:3
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

    % fixa escala desde o começo
    set(ax, 'XLim', [X1(1) X1(end)], 'XLimMode','manual');
    set(ax, 'YLim', YLimFixed,       'YLimMode','manual');

    % acha linhas por DisplayName
    plines = findall(ax, "Type","line");
    dn = string(get(plines,"DisplayName"));

    hs = gobjects(1,3);
    for j = 1:3
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
    for j = 1:3
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
