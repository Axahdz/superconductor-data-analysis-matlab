% Superconductor Data Analysis Script
%
% This script processes experimental superconducting data files.
%
% Requirements:
% The data files (.Cn1 and .Cn2) must be located in the same
% folder as this script.
%
% Example data can be found in the folder:
% example_dataset
%
% Author: Javier Axayácatl Melchor Hernández
%% ================================================================
%  Doble Ajuste (Lorentziano y Gaussiano) — Arreglo Lineal
%  Selección: 'lorentz' | 'gauss' | 'both'
% ================================================================
clear; clc; close all;

%% ---------------- CONFIGURACIÓN DEL EXPERIMENTO ------------------------
% Ruta de los archivos Run0_*.Cn1 (Va) y Run0_*.Cn2 (Vb)
ruta = fileparts(mfilename('fullpath'));

% Experimentos (lineal): usa los mismos IDs en ambos modelos
expIDs  = {'01968','01966','01964','01962','01970','01972'};

% Letras para etiquetar subplots
letters = 'abcdefghijklmnopqrstuvwxyz';

% Qué correr: 'lorentz' | 'gauss' | 'both'
model_to_run = 'both';

% Colores consistentes con tus scripts
colDV  = [0.2 0.6 0.2];   % puntos/derivada
colG1  = [0.1 0.6 1.0];   % azul (curva 1)
colG2  = [1.0 0.6 0.0];   % naranja (curva 2)
colSum = [0.8 0.2 0.2];   % rojo  (suma)

% Carpetas de salida
outBase = fileparts(mfilename('fullpath'));
outDirLorentz = fullfile(outBase,'2_lorentzianas');
outDirGauss   = fullfile(outBase,'2_gaussianas');
if ~exist(outDirLorentz,'dir'), mkdir(outDirLorentz); end
if ~exist(outDirGauss,'dir'),   mkdir(outDirGauss);   end

%% ------------ PASO 1: extraer y ordenar por corriente ------------------
I_vect = zeros(1,numel(expIDs));
for n = 1:numel(expIDs)
    fid = fopen(fullfile(ruta,['Run0_',expIDs{n},'.Cn1']),'r');
    firstLine = fgetl(fid); fclose(fid);
    vals = textscan(firstLine,'%f','Delimiter',',');
    I_vect(n) = vals{1}(3)*1e3;  % A -> mA
end
[I_vect, idx] = sort(I_vect,'ascend');
expIDs = expIDs(idx);

%% ------------------- EJECUCIÓN SEGÚN SELECCIÓN -------------------------
switch lower(model_to_run)
    case 'lorentz'
        run_lorentz(ruta, expIDs, I_vect, letters, colDV, colG1, colG2, colSum, outDirLorentz);
    case 'gauss'
        run_gauss(ruta, expIDs, I_vect, letters, colDV, colG1, colG2, colSum, outDirGauss);
    case 'both'
        run_lorentz(ruta, expIDs, I_vect, letters, colDV, colG1, colG2, colSum, outDirLorentz);
        run_gauss(ruta, expIDs, I_vect, letters, colDV, colG1, colG2, colSum, outDirGauss);
    otherwise
        error('model_to_run debe ser ''lorentz'', ''gauss'' o ''both''.');
end

disp('--- Proceso terminado ---');

%% ======================= FUNCIONES LOCALES =============================
function run_lorentz(ruta, expIDs, I_vect, letters, colDV, colG1, colG2, colSum, outDir)
    % ----------------- Panel figura -----------------
    nRows = ceil(numel(expIDs)/2);
    fig = figure('Units','normalized','Position',[0.20 0.05 0.70 0.18+0.25*nRows],'Color','w','Name','Panel doble Lorentziana (lineal)');
    tl = tiledlayout(nRows,2,'Padding','compact','TileSpacing','compact','OuterPosition',[0.02 0.05 0.85 0.9]);
    xlabel(tl,'Temperatura (K)','FontSize',23,'FontWeight','bold');
    ylabel(tl,'dV_{prom}/dT (\muV/K)','FontSize',28,'FontWeight','bold','Color',colDV);

    % ----------------- Tabla resultados -----------------
    resultados = table('Size',[numel(expIDs)*2 14], ...
        'VariableTypes',{'string','double','double','double','double','double','double','double','double','double','string','double','double','string'}, ...
        'VariableNames',{'Experimento','Corriente_mA','Amplitud_A','Centro_x0','Gamma','Ancho_gamma','Area_lorentz','Varianza','StdDev','Max_dVdT','Componente','CentroTotal_x0','Rsquare','AjusteOK'});

    row = 1;

    for k = 1:numel(expIDs)
        id = expIDs{k};
        [Temp, dVdT] = cargar_datos_y_derivada(ruta, id);

        % Modelo doble Lorentziana
        modelo = fittype('A1*(0.5*s1)/((x-x1)^2 + (0.5*s1)^2) + A2*(0.5*s2)/((x-x2)^2 + (0.5*s2)^2)', ...
            'independent','x','coefficients',{'A1','x1','s1','A2','x2','s2'});

        try
            [Amax, imax] = max(dVdT); x_est = Temp(imax);
            [ajuste, gof] = fit(Temp, dVdT, modelo, ...
                'StartPoint',[Amax*0.2, x_est-1, 0.3, Amax*0.8, x_est, 0.5], ...
                'Lower',[0, x_est-1.5, 0.05, 0, x_est-0.3, 0.05], ...
                'Upper',[Inf, x_est+1.0, 2.0, Inf, x_est+0.3, 2.0]);

            fitX = linspace(min(Temp), max(Temp), 800);
            g1 = ajuste.A1*(0.5*ajuste.s1)./((fitX-ajuste.x1).^2 + (0.5*ajuste.s1)^2);
            g2 = ajuste.A2*(0.5*ajuste.s2)./((fitX-ajuste.x2).^2 + (0.5*ajuste.s2)^2);
            gT = g1 + g2;
            [~, idxMax] = max(gT); centroTotal = fitX(idxMax);

            area1 = pi * ajuste.A1 * ajuste.s1;
            area2 = pi * ajuste.A2 * ajuste.s2;
            var1 = (ajuste.s1/2)^2 * (pi/2);
            var2 = (ajuste.s2/2)^2 * (pi/2);
            std1 = sqrt(var1); std2 = sqrt(var2);
            max1 = max(g1);   max2 = max(g2);

            resultados(row,:)   = {id, I_vect(k), ajuste.A1, ajuste.x1, ajuste.s1, ajuste.s1, area1, var1, std1, max1, 'Lorentz1', centroTotal, gof.rsquare, 'Sí'};
            resultados(row+1,:) = {id, I_vect(k), ajuste.A2, ajuste.x2, ajuste.s2, ajuste.s2, area2, var2, std2, max2, 'Lorentz2', centroTotal, gof.rsquare, 'Sí'};

            % ---- Gráfica por experimento (panel) ----
            ax = nexttile(tl); hold(ax,'on'); set(ax,'Color','w');
            plot(ax, Temp, dVdT, 'o', 'MarkerFaceColor', colDV, 'MarkerEdgeColor', colDV*0.8);
            plot(ax, fitX, g1, '--','Color',colG1,'LineWidth',2.0);
            plot(ax, fitX, g2, '--','Color',colG2,'LineWidth',2.0);
            plot(ax, fitX, gT,  '-','Color',colSum,'LineWidth',1.6);
            if k==2
                legend(ax, {'Datos','Lorentz 1','Lorentz 2','Suma'}, 'FontSize',10,'Box','on','Location','northeast');
            end
            estilizar_axes_lineal(ax, Temp, dVdT, letters, k);

            title(ax, ['I = ',num2str(I_vect(k),'%.1f'),' mA (Exp ',id,')'], 'FontSize',16,'FontWeight','bold');
            row = row + 2;

        catch ME
            warning(['Fallo ajuste Lorentz en ', id, ': ', ME.message]);
            resultados(row,:)   = {id, I_vect(k), NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'Lorentz1', NaN, NaN, 'No'};
            resultados(row+1,:) = {id, I_vect(k), NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'Lorentz2', NaN, NaN, 'No'};
            row = row + 2;
        end
    end

    % ------------------- Exportación Excel/PNG -------------------
    resultadosFile = fullfile(outDir,'resultados_doble_lorentziana_lineal.xlsx');

    TablaAzul    = resultados(strcmp(resultados.Componente, 'Lorentz1'), :);
    TablaNaranja = resultados(strcmp(resultados.Componente, 'Lorentz2'), :);

    nA = height(TablaAzul); nN = height(TablaNaranja);
    Corriente_fmt_A = strings(nA,1); Amplitud_fmt_A = strings(nA,1); Centro_fmt_A = strings(nA,1); Gamma_fmt_A = strings(nA,1); Tc_fmt_A = strings(nA,1);
    Corriente_fmt_N = strings(nN,1); Amplitud_fmt_N = strings(nN,1); Centro_fmt_N = strings(nN,1); Gamma_fmt_N = strings(nN,1); Tc_fmt_N = strings(nN,1);

    for i = 1:nA
        I = TablaAzul.Corriente_mA(i); A = TablaAzul.Amplitud_A(i); x0 = TablaAzul.Centro_x0(i); G = TablaAzul.Gamma(i);
        Corriente_fmt_A(i) = sprintf('%.1f ± 0.0001', I);
        Amplitud_fmt_A(i)  = sprintf('%.3f ± 0.001', A);
        Centro_fmt_A(i)    = sprintf('%.3f', x0);
        Gamma_fmt_A(i)     = sprintf('%.2f', G);
        Tc_fmt_A(i)        = sprintf('%.1f ± %.1f', round(x0,1), round(G/2,1));
    end
    for i = 1:nN
        I = TablaNaranja.Corriente_mA(i); A = TablaNaranja.Amplitud_A(i); x0 = TablaNaranja.Centro_x0(i); G = TablaNaranja.Gamma(i);
        Corriente_fmt_N(i) = sprintf('%.1f ± 0.0001', I);
        Amplitud_fmt_N(i)  = sprintf('%.3f ± 0.001', A);
        Centro_fmt_N(i)    = sprintf('%.3f', x0);
        Gamma_fmt_N(i)     = sprintf('%.2f', G);
        Tc_fmt_N(i)        = sprintf('%.1f ± %.1f', round(x0,1), round(G/2,1));
    end

    TablaAzul_fmt = table(Corriente_fmt_A, TablaAzul.Experimento, Amplitud_fmt_A, Centro_fmt_A, Gamma_fmt_A, Tc_fmt_A, ...
        TablaAzul.Area_lorentz, TablaAzul.Varianza, TablaAzul.StdDev, ...
        'VariableNames', {'Corriente','Experimento','Amplitud','Centro','Gamma','Tc','Area','Varianza','StdDev'});

    TablaNaranja_fmt = table(Corriente_fmt_N, TablaNaranja.Experimento, Amplitud_fmt_N, Centro_fmt_N, Gamma_fmt_N, Tc_fmt_N, ...
        TablaNaranja.Area_lorentz, TablaNaranja.Varianza, TablaNaranja.StdDev, ...
        'VariableNames', {'Corriente','Experimento','Amplitud','Centro','Gamma','Tc','Area','Varianza','StdDev'});

    experimentos = unique(resultados.Experimento);
    tablaTotal = table('Size',[numel(experimentos) 6], ...
        'VariableTypes',{'string','string','double','double','string','double'}, ...
        'VariableNames',{'Corriente','Experimento','Area_total','CentroMax_total','T_c','Rsquare'});

    for i = 1:numel(experimentos)
        id = experimentos{i};
        datos = resultados(strcmp(resultados.Experimento,id),:);
        I_fmt = sprintf('%.1f ± 0.0001', datos.Corriente_mA(1));
        if all(strcmp(datos.AjusteOK, 'Sí'))
            centro = datos.CentroTotal_x0(1);
            Tc_str = sprintf('%.1f ± 0.5', round(centro,1));
            tablaTotal(i,:) = {I_fmt, id, sum(datos.Area_lorentz), centro, Tc_str, datos.Rsquare(1)};
        else
            tablaTotal(i,:) = {I_fmt, id, NaN, NaN, "", NaN};
        end
    end

    writetable(TablaAzul_fmt, resultadosFile, 'Sheet','Curva azul (Lorentz1)');
    writetable(TablaNaranja_fmt, resultadosFile, 'Sheet','Curva naranja (Lorentz2)');
    writetable(tablaTotal, resultadosFile, 'Sheet','Curva roja total (Suma)');
    disp(['[Lorentz] Excel: ', resultadosFile]);

    outFile = fullfile(outDir, ['panel_doble_lorentziana_lineal',expIDs{1},'-',expIDs{end},'.png']);
    exportgraphics(fig, outFile, 'Resolution', 300, 'BackgroundColor', 'white');
    disp(['[Lorentz] Panel PNG: ', outFile]);

    % ---- Gráficas individuales (azul, naranja, total) ----
    figL1 = figure('Name','Corriente vs Centro - Lorentz1','Color','w','Units','normalized','Position',[0.3 0.3 0.4 0.4]);
    plot(TablaAzul.Corriente_mA, TablaAzul.Centro_x0, 'o-', 'Color', colG1, 'MarkerFaceColor', colG1, 'MarkerEdgeColor', 'k', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Corriente (mA)', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Centro Lorentz1 (K)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Curva azul (Lorentz1)', 'FontSize', 16, 'FontWeight', 'bold'); grid on; box on;
    exportgraphics(figL1, fullfile(outDir, 'corriente_vs_centro_Lorentz1.png'), 'Resolution', 300);

    figL2 = figure('Name','Corriente vs Centro - Lorentz2','Color','w','Units','normalized','Position',[0.35 0.35 0.4 0.4]);
    plot(TablaNaranja.Corriente_mA, TablaNaranja.Centro_x0, 'o-', 'Color', colG2, 'MarkerFaceColor', colG2, 'MarkerEdgeColor', 'k', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Corriente (mA)', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Centro Lorentz2 (K)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Curva naranja (Lorentz2)', 'FontSize', 16, 'FontWeight', 'bold'); grid on; box on;
    exportgraphics(figL2, fullfile(outDir, 'corriente_vs_centro_Lorentz2.png'), 'Resolution', 300);

    Corr_mA_total = str2double(extractBefore(tablaTotal.Corriente, ' ±'));
    figLT = figure('Name','Corriente vs Centro - Lorentz total','Color','w','Units','normalized','Position',[0.4 0.4 0.4 0.4]);
    plot(Corr_mA_total, tablaTotal.CentroMax_total, 'o-', 'Color', colSum, 'MarkerFaceColor', colSum, 'MarkerEdgeColor', 'k', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Corriente (mA)', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Centro máx total (K)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Curva roja total (Lorentz1 + Lorentz2)', 'FontSize', 16, 'FontWeight', 'bold'); grid on; box on;
    exportgraphics(figLT, fullfile(outDir, 'corriente_vs_centro_Total_Lorentz.png'), 'Resolution', 300);

    % ---- Tabla Tc Azul vs Naranja + gráfica comparativa con errores ----
    Tc_azul    = strings(nA,1);
    Tc_naranja = strings(nN,1);
    for i = 1:nA
        x0a = TablaAzul.Centro_x0(i); Ga = TablaAzul.Gamma(i);
        Tc_azul(i) = sprintf('%.1f ± %.1f', round(x0a,1), round(Ga/2,1));
    end
    for i = 1:nN
        x0n = TablaNaranja.Centro_x0(i); Gn = TablaNaranja.Gamma(i);
        Tc_naranja(i) = sprintf('%.1f ± %.1f', round(x0n,1), round(Gn/2,1));
    end
    tabla_Tc_comparada = table( ...
        Corriente_fmt_A, TablaAzul.Experimento, Tc_azul, Tc_naranja, ...
        'VariableNames', {'Corriente','Experimento','T_c_Curva_azul','T_c_Curva_naranja'});

    writetable(tabla_Tc_comparada, resultadosFile, 'Sheet','T_c Azul y Naranja');

    figComp = figure('Name','Comparación Tc (Lorentz) con errores','Color','w','Units','normalized','Position',[0.35 0.35 0.45 0.45]);
    hold on;
    TcAzul     = TablaAzul.Centro_x0;  ErrAzul    = TablaAzul.Gamma / 2;
    TcNaranja  = TablaNaranja.Centro_x0; ErrNaranja = TablaNaranja.Gamma / 2;
    errorbar(TablaNaranja.Corriente_mA, TcNaranja, ErrNaranja, '-o', 'Color', colG2, 'MarkerFaceColor', colG2, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'CapSize', 10);
    errorbar(TablaAzul.Corriente_mA,   TcAzul,    ErrAzul,    '-s', 'Color', colG1, 'MarkerFaceColor', colG1, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'CapSize', 10);
    xlabel('Corriente (mA)', 'FontSize', 22, 'FontWeight', 'bold');
    ylabel('Temperatura Crítica (K)', 'FontSize', 22, 'FontWeight', 'bold');
    legend({'T_c Principal','T_c Secundaria'}, 'FontSize', 12, 'Location','northeast');
    grid on; box on; set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    xlim([0.4 1.6]); xticks(0.4:0.1:1.6);
    ylim([86 91]);   yticks(86:1:91);
    exportgraphics(figComp, fullfile(outDir, 'Tc_comparacion_con_errores.png'), 'Resolution', 300, 'BackgroundColor', 'white');
end

function run_gauss(ruta, expIDs, I_vect, letters, colDV, colG1, colG2, colSum, outDir)
    % ----------------- Panel figura -----------------
    nRows = ceil(numel(expIDs)/2);
    fig = figure('Units','normalized','Position',[0.20 0.05 0.70 0.18+0.25*nRows],'Color','w','Name','Panel doble Gaussiana (lineal)');
    tl = tiledlayout(nRows,2,'Padding','compact','TileSpacing','compact','OuterPosition',[0.02 0.05 0.85 0.9]);
    xlabel(tl,'Temperatura (K)','FontSize',23,'FontWeight','bold');
    ylabel(tl,'dV_{prom}/dT (\muV/K)','FontSize',28,'FontWeight','bold','Color',colDV);

    % ----------------- Tabla resultados -----------------
    resultados = table('Size',[numel(expIDs)*2 14], ...
        'VariableTypes',{'string','double','double','double','double','double','double','double','double','double','string','double','double','string'}, ...
        'VariableNames',{'Experimento','Corriente_mA','Amplitud_A','Centro_x0','Ancho_sigma','FWHM','Area_gaussiana','Varianza','StdDev','Max_dVdT','Componente','CentroTotal_x0','Rsquare','AjusteOK'});

    row = 1;

    for k = 1:numel(expIDs)
        id = expIDs{k};
        [Temp, dVdT] = cargar_datos_y_derivada(ruta, id);

        % Modelo doble Gaussiana
        modelo = fittype('A1*exp(-((x-x1)^2)/(2*s1^2)) + A2*exp(-((x-x2)^2)/(2*s2^2))', ...
            'independent','x', 'coefficients',{'A1','x1','s1','A2','x2','s2'});

        try
            [Amax, imax] = max(dVdT); x_est = Temp(imax);
            x1_0 = x_est - 1.0;  x2_0 = x_est;
            A1_0 = Amax * 0.2;   A2_0 = Amax * 0.8;
            s1_0 = 0.3;          s2_0 = 0.5;

            [ajuste, gof] = fit(Temp, dVdT, modelo, ...
                'StartPoint', [A1_0, x1_0, s1_0, A2_0, x2_0, s2_0], ...
                'Lower',      [0, x1_0 - 0.5, 0.05, 0, x2_0 - 0.3, 0.05], ...
                'Upper',      [Inf, x1_0 + 1.0, 2.0,  Inf, x2_0 + 0.3, 2.0]);

            fitX = linspace(min(Temp), max(Temp), 800);
            g1 = ajuste.A1 * exp(-((fitX - ajuste.x1).^2) / (2 * ajuste.s1^2));
            g2 = ajuste.A2 * exp(-((fitX - ajuste.x2).^2) / (2 * ajuste.s2^2));
            gT = g1 + g2;
            [~, idxMax] = max(gT); centroTotal = fitX(idxMax);

            area1 = ajuste.A1 * ajuste.s1 * sqrt(2*pi);
            area2 = ajuste.A2 * ajuste.s2 * sqrt(2*pi);
            var1 = ajuste.s1^2;  std1 = ajuste.s1;
            var2 = ajuste.s2^2;  std2 = ajuste.s2;
            max1 = max(g1);      max2 = max(g2);

            resultados(row,:)   = {id, I_vect(k), ajuste.A1, ajuste.x1, ajuste.s1, 2.3548*ajuste.s1, area1, var1, std1, max1, 'Gauss1', centroTotal, gof.rsquare, 'Sí'};
            resultados(row+1,:) = {id, I_vect(k), ajuste.A2, ajuste.x2, ajuste.s2, 2.3548*ajuste.s2, area2, var2, std2, max2, 'Gauss2', centroTotal, gof.rsquare, 'Sí'};

            % ---- Gráfica por experimento (panel) ----
            ax = nexttile(tl); set(ax,'Color','w'); hold(ax,'on');
            plot(ax, Temp, dVdT, 'o', 'MarkerFaceColor', colDV, 'MarkerEdgeColor', colDV*0.8, 'MarkerSize', 5, 'LineWidth', 1.2);
            plot(ax, fitX, g1, '--','Color',colG1,'LineWidth',2.0);
            plot(ax, fitX, g2, '--','Color',colG2,'LineWidth',2.0);
            plot(ax, fitX, gT,  '-','Color',colSum,'LineWidth',1.6);
            if k==2
                legend(ax, {'Datos','Gauss 1','Gauss 2','Suma'}, 'FontSize',10,'Box','on','Location','northeast');
            end
            estilizar_axes_lineal(ax, Temp, dVdT, letters, k);
            title(ax, ['I = ',num2str(I_vect(k),'%.1f'),' mA  (Exp ',id,')'],'FontSize',16,'FontWeight','bold');

            row = row + 2;

        catch ME
            warning(['Fallo ajuste Gauss en ', id, ': ', ME.message]);
            resultados(row,:)   = {id, I_vect(k), NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'Gauss1', NaN, NaN, 'No'};
            resultados(row+1,:) = {id, I_vect(k), NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'Gauss2', NaN, NaN, 'No'};
            row = row + 2;
        end
    end

    % ------------------- Exportación Excel/PNG -------------------
    resultadosFile = fullfile(outDir,'resultados_doble_gaussiana_lineal.xlsx');

    TablaAzul    = resultados(strcmp(resultados.Componente, 'Gauss1'), :);
    TablaNaranja = resultados(strcmp(resultados.Componente, 'Gauss2'), :);

    nA = height(TablaAzul); nN = height(TablaNaranja);
    Corriente_fmt_A = strings(nA,1); Amplitud_fmt_A = strings(nA,1); Centro_fmt_A = strings(nA,1); FWHM_fmt_A = strings(nA,1); Tc_fmt_A = strings(nA,1);
    Corriente_fmt_N = strings(nN,1); Amplitud_fmt_N = strings(nN,1); Centro_fmt_N = strings(nN,1); FWHM_fmt_N = strings(nN,1); Tc_fmt_N = strings(nN,1);

    for i = 1:nA
        I = TablaAzul.Corriente_mA(i); A = TablaAzul.Amplitud_A(i); x0 = TablaAzul.Centro_x0(i); FWHM = TablaAzul.FWHM(i);
        Corriente_fmt_A(i) = sprintf('%.1f ± 0.0001', I);
        Amplitud_fmt_A(i)  = sprintf('%.3f ± 0.001', A);
        Centro_fmt_A(i)    = sprintf('%.3f', x0);
        FWHM_fmt_A(i)      = sprintf('%.2f', FWHM);
        Tc_fmt_A(i)        = sprintf('%.1f ± %.1f', round(x0,1), round(0.5*FWHM,1));
    end
    for i = 1:nN
        I = TablaNaranja.Corriente_mA(i); A = TablaNaranja.Amplitud_A(i); x0 = TablaNaranja.Centro_x0(i); FWHM = TablaNaranja.FWHM(i);
        Corriente_fmt_N(i) = sprintf('%.1f ± 0.0001', I);
        Amplitud_fmt_N(i)  = sprintf('%.3f ± 0.001', A);
        Centro_fmt_N(i)    = sprintf('%.3f', x0);
        FWHM_fmt_N(i)      = sprintf('%.2f', FWHM);
        Tc_fmt_N(i)        = sprintf('%.1f ± %.1f', round(x0,1), round(0.5*FWHM,1));
    end

    TablaAzul_fmt = table(Corriente_fmt_A, TablaAzul.Experimento, Amplitud_fmt_A, Centro_fmt_A, FWHM_fmt_A, Tc_fmt_A, ...
        TablaAzul.Ancho_sigma, TablaAzul.Area_gaussiana, TablaAzul.Varianza, TablaAzul.StdDev, ...
        'VariableNames', {'Corriente','Experimento','Amplitud','Centro','FWHM','Tc','Sigma','Area','Varianza','StdDev'});

    TablaNaranja_fmt = table(Corriente_fmt_N, TablaNaranja.Experimento, Amplitud_fmt_N, Centro_fmt_N, FWHM_fmt_N, Tc_fmt_N, ...
        TablaNaranja.Ancho_sigma, TablaNaranja.Area_gaussiana, TablaNaranja.Varianza, TablaNaranja.StdDev, ...
        'VariableNames', {'Corriente','Experimento','Amplitud','Centro','FWHM','Tc','Sigma','Area','Varianza','StdDev'});

    experimentos = unique(resultados.Experimento);
    tablaTotal = table('Size',[numel(experimentos) 6], ...
        'VariableTypes',{'string','string','double','double','string','double'}, ...
        'VariableNames',{'Corriente','Experimento','Area_total','CentroMax_total','T_c','Rsquare'});

    for i = 1:numel(experimentos)
        id = experimentos{i};
        datos = resultados(strcmp(resultados.Experimento,id),:);
        I_fmt = sprintf('%.1f ± 0.0001', datos.Corriente_mA(1));
        if all(strcmp(datos.AjusteOK, 'Sí'))
            centro = datos.CentroTotal_x0(1);
            Tc_str = sprintf('%.1f ± 0.5', round(centro,1));
            tablaTotal(i,:) = {I_fmt, id, sum(datos.Area_gaussiana), centro, Tc_str, datos.Rsquare(1)};
        else
            tablaTotal(i,:) = {I_fmt, id, NaN, NaN, "", NaN};
        end
    end

    writetable(TablaAzul_fmt, resultadosFile, 'Sheet','Curva azul (Gauss1)');
    writetable(TablaNaranja_fmt, resultadosFile, 'Sheet','Curva naranja (Gauss2)');
    writetable(tablaTotal, resultadosFile, 'Sheet','Curva roja total (Suma)');
    disp(['[Gauss] Excel: ', resultadosFile]);

    outFile = fullfile(outDir, ['panel_doble_gaussiana_lineal',expIDs{1},'-',expIDs{end},'.png']);
    exportgraphics(fig, outFile, 'Resolution', 300, 'BackgroundColor', 'white');
    disp(['[Gauss] Panel PNG: ', outFile]);

    % ---- Gráficas individuales (azul, naranja, total) ----
    figG1 = figure('Name','Corriente vs Centro - Gauss1','Color','w','Units','normalized','Position',[0.3 0.3 0.4 0.4]);
    plot(TablaAzul.Corriente_mA, TablaAzul.Centro_x0, 'o-', 'Color', colG1, 'MarkerFaceColor', colG1, 'MarkerEdgeColor', 'k', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Corriente (mA)', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Centro Gauss1 (K)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Curva azul (Gauss1)', 'FontSize', 16, 'FontWeight', 'bold'); grid on; box on;
    exportgraphics(figG1, fullfile(outDir, 'corriente_vs_centro_Gauss1.png'), 'Resolution', 300);

    figG2 = figure('Name','Corriente vs Centro - Gauss2','Color','w','Units','normalized','Position',[0.35 0.35 0.4 0.4]);
    plot(TablaNaranja.Corriente_mA, TablaNaranja.Centro_x0, 'o-', 'Color', colG2, 'MarkerFaceColor', colG2, 'MarkerEdgeColor', 'k', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Corriente (mA)', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Centro Gauss2 (K)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Curva naranja (Gauss2)', 'FontSize', 16, 'FontWeight', 'bold'); grid on; box on;
    exportgraphics(figG2, fullfile(outDir, 'corriente_vs_centro_Gauss2.png'), 'Resolution', 300);

    Corr_mA_total = str2double(extractBefore(tablaTotal.Corriente, ' ±'));
    figGT = figure('Name','Corriente vs Centro - Gauss total','Color','w','Units','normalized','Position',[0.4 0.4 0.4 0.4]);
    plot(Corr_mA_total, tablaTotal.CentroMax_total, 'o-', 'Color', colSum, 'MarkerFaceColor', colSum, 'MarkerEdgeColor', 'k', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Corriente (mA)', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Centro máx total (K)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Curva roja total (Gauss1 + Gauss2)', 'FontSize', 16, 'FontWeight', 'bold'); grid on; box on;
    exportgraphics(figGT, fullfile(outDir, 'corriente_vs_centro_Total_Gauss.png'), 'Resolution', 300);
end

function [Temp, dVdT] = cargar_datos_y_derivada(ruta, id)
    % Lee Va/Vb, calcula Vprom y dV/dT, y recorta al rango [80, 100] K
    fVA = fullfile(ruta,['Run0_',id,'.Cn1']);
    fVB = fullfile(ruta,['Run0_',id,'.Cn2']);
    opts = {'Delimiter',',','ReadVariableNames',false,'FileType','text'};
    Tva = readtable(fVA,opts{:});
    Tvb = readtable(fVB,opts{:});
    Tva.Properties.VariableNames = {'Tiempo','T','I','V','Err','Cfg'};
    Tvb.Properties.VariableNames = {'Tiempo','T','I','V','Err','Cfg'};

    Vprom = 0.5 * (Tva.V - Tvb.V) * 1e6;    % µV
    Temp  = 0.5*(Tva.T + Tvb.T);            % K

    sel   = Temp >= 80 & Temp <= 100;
    Temp  = Temp(sel);
    Vprom = Vprom(sel);
    if numel(Temp) < 10
        warning(['Muy pocos puntos en ', id]);
    end

    dVdT = gradient(Vprom) ./ gradient(Temp);
end

function estilizar_axes_lineal(ax, Temp, dVdT, letters, k)
    % Estilo de ejes consistente con tus scripts (rango X, Y y marca de letra)
    grid(ax,'on'); grid(ax,'minor');
    xlim(ax,[82 98]); xticks(ax,82:2:98);

    % Reglas especiales según letra (como tu script lorentziano)
    switch letters(k)
        case 'a'
            ylim(ax, [0 5]);   yticks(ax, 0:1:5);
        case 'e'
            ylim(ax, [0 12]);  yticks(ax, 0:2:12);
        otherwise
            ymax = ceil(max(dVdT));
            if ymax == 7, ymax = 8; end
            if ymax < 2, ymax = 2; end
            ylim(ax, [0 ymax]);
            if ymax < 6
                yticks(ax,0:1:ymax);
            elseif mod(ymax,2)==0
                yticks(ax,0:2:ymax);
            else
                yticks(ax,0:3:ymax);
            end
    end

    text(ax,0.08,0.85,['(' letters(k) ')'],'Units','normalized','FontWeight','bold','FontSize',24);
    set(ax,'FontSize',12,'FontWeight','bold','LineWidth',1);
end
