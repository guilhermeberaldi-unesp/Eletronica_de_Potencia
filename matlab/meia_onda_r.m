%% Fonte senoidal v_s(t) + retificador de meia-onda (1 diodo) + carga R


clear; clc; close all;

% ---- Parâmetros ----
f   = 60;          % Hz
Vm  = 170;         % Vp (pico)
R   = 10;          % ohms
Vd  = 0.7;         % V (queda do diodo; use 0 para ideal)
ncy = 3;           % número de ciclos
fs  = 200e3;       % Hz

% ---- Tempo ----
T = 1/f;
t = 0:1/fs:(ncy*T);

% ---- Sinais ----
v_s = Vm*sin(2*pi*f*t);     % v_s(t)
v_0 = max(v_s - Vd, 0);     % v_0(t)
i_0 = v_0 / R;              % i_0(t)

% Tensão sobre o diodo v_D(t) (anodo->catodo):
% quando conduz: ~Vd; quando bloqueia: v_s(t)
v_D = v_s - v_0;

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

% ---- Estilo de linha ----
LW = 3.0;

% ---- Eixo x em radianos ----
theta = 2*pi*f*t;   % rad

figure('Color','w');
plot(theta, v_s, 'Color', plotColors(1,:), 'LineWidth', LW); hold on;
plot(theta, v_0,'Color', plotColors(2,:), 'LineWidth', LW);
plot(theta, i_0, 'Color', plotColors(3,:), 'LineWidth', LW);
plot(theta, v_D,'--', 'Color', plotColors(4,:), 'LineWidth', LW);  % próxima cor (red)
grid on;

xlabel('\theta (rad)');
ylabel('Tensão (V) / Corrente (A)');
title('Retificador monofásico de meia onda com carga resistiva');

% Ticks em múltiplos de pi
xlim([0, ncy*2*pi]);
xticks(0:pi:(ncy*2*pi));
xticklabels(compose('%g\\pi', (0:1:(ncy*2))));
xticklabels(strrep(xticklabels,'0\pi','0'));

% Legenda 
lgd = legend('v_s(t)','v_0(t)','i_0(t)','v_D(t)','Location','best');
lgd.FontSize   = 12;
lgd.FontWeight = 'bold';
lgd.Color      = 'none';

% Fontes dos eixos em negrito
ax = gca;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontWeight = 'bold';

%% =========================
%  ADENDO: 4 MP4s (um por linha), usando a .fig como TEMPLATE
%  V1: anima v_s
%  V2: v_s fixa + anima v_0
%  V3: v_s,v_0 fixas + anima i_0
%  V4: v_s,v_0,i_0 fixas + anima v_D
% =========================

figTemplate = "meia onda r.fig";          % no mesmo folder

% Dados
X1 = theta(:);
Y = [v_s(:), v_0(:), i_0(:), v_D(:)];
N = numel(X1);

if ~isfile(figTemplate)
    error("Não achei '%s' em %s", figTemplate, pwd);
end

% Config do vídeo
fps      = 30;
durEach  = 10;                    % duração de CADA vídeo (s)
nFrames  = max(60, round(durEach*fps));
resDPI   = 300;                   % pode pôr 600 (mais lento)

% Helper
reveal = @(y,idx) [y(1:idx); nan(N-idx,1)];

% Nomes e DisplayNames na ordem
names = ["v_s(t)","v_0(t)","i_0(t)","v_D(t)"];
outFiles = ["meia_onda_v1_vs.mp4", ...
            "meia_onda_v2_v0.mp4", ...
            "meia_onda_v3_i0.mp4", ...
            "meia_onda_v4_vD.mp4"];

for vid = 1:4
    % Abre template novo para cada vídeo (evita "herdar" estado)
    fig = openfig(figTemplate, "new", "visible");
    drawnow;

    % trava tamanho (opcional: defina Position se quiser 1080p/4K)
    set(fig,'Units','pixels');
    set(fig,'Resize','off');
    drawnow;

    % Acha axes principal
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

    % Acha linhas por DisplayName
    plines = findall(ax, "Type","line");
    getByName = @(nm) plines(strcmp(string(get(plines,"DisplayName")), nm));

    h1 = getByName("v_s(t)");
    h2 = getByName("v_0(t)");
    h3 = getByName("i_0(t)");
    h4 = getByName("v_D(t)");

    if isempty(h1) || isempty(h2) || isempty(h3) || isempty(h4)
        error("Não achei as 4 curvas pelo DisplayName dentro da FIG.");
    end

    hs = [h1 h2 h3 h4];

    % Garante XData certo
    set(hs, "XData", X1);

    % Estado inicial para ESTE vídeo:
    % - curvas anteriores (1..vid-1) completas
    % - curva vid vazia (vai animar)
    % - curvas futuras vazias
    for j = 1:4
        if j < vid
            set(hs(j), "YData", Y(:,j));
        else
            set(hs(j), "YData", nan(N,1));
        end
    end
    drawnow;

    % Prepara VideoWriter
    vw = VideoWriter(outFiles(vid), "MPEG-4");
    vw.FrameRate = fps;
    vw.Quality   = 100;
    open(vw);
    cleanupObj = onCleanup(@() close(vw));

    % Temp png único por vídeo (evita conflitos)
    tmpPng = fullfile(tempdir, "frame_tmp_" + string(feature("getpid")) + "_v" + vid + ".png");

    % Frame "teste" para fixar dimensão
    exportgraphics(fig, tmpPng, 'Resolution', resDPI);
    im0 = imread(tmpPng);
    H = size(im0,1); W = size(im0,2);

    % Anima SOMENTE a curva vid
    for k = 1:nFrames
        frac = (k-1)/(nFrames-1);
        idx  = max(2, round(frac*(N-1)) + 1);

        set(hs(vid), "YData", reveal(Y(:,vid), idx));
        drawnow;

        exportgraphics(fig, tmpPng, 'Resolution', resDPI);
        im = imread(tmpPng);

        % segurança: garante tamanho fixo
        if size(im,1) ~= H || size(im,2) ~= W
            im = im(1:min(end,H), 1:min(end,W), :);
        end

        writeVideo(vw, im);
    end

    % Fecha e limpa
    if exist(tmpPng,'file'), delete(tmpPng); end
    clear cleanupObj;   % close(vw)

    close(fig);

    disp("Salvou: " + fullfile(pwd, outFiles(vid)));
end