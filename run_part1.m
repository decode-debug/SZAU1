%% SZAU, zadanie 1, czesc 1: zbiorniki i konwencjonalny DMC
% Uruchom z folderu projektu: run_part1
% Wymaga Control System Toolbox (ss, c2d, step).
clear; close all; clc;
set(groot,'defaultFigureColor','w','defaultAxesColor','w', ...
    'defaultAxesXColor','k','defaultAxesYColor','k', ...
    'defaultTextColor','k','defaultLegendColor','w', ...
    'defaultLegendTextColor','k');
%% Parametry i punkt pracy
P.A1 = 560;       % pole przekroju zbiornika 1 [cm^2]
P.C2 = 0.6;       % V2 = C2*h2^3
P.a1 = 29;        % wspolczynnik odplywu ze zbiornika 1
P.a2 = 33;        % wspolczynnik odplywu ze zbiornika 2
P.u0 = 100;       % F1in w punkcie pracy [cm^3/s]
P.d0 = 20;        % FD w punkcie pracy [cm^3/s]
P.tau = 200;      % opoznienie sterowania [s]
P.Ts = 20;        % okres probkowania [s]
P.delay = P.tau/P.Ts;
P.q0 = P.u0 + P.d0;
P.h10 = (P.q0/P.a1)^2;
P.h20 = (P.q0/P.a2)^2;
P.v10 = P.A1*P.h10;
P.v20 = P.C2*P.h20^3;

%% Model liniowy w postaci state-space
P.k1 = P.a1/(2*sqrt(P.h10));
P.k2 = P.a2/(2*sqrt(P.h20));
P.cap2 = 3*P.C2*P.h20^2;
P.Ac = [-P.k1/P.A1, 0; P.k1/P.cap2, -P.k2/P.cap2];
P.Bc = [1/P.A1, 1/P.A1; 0, 0]; % wejscia: Delta F1in, Delta FD
P.C = [0, 1];
P.sys = ss(P.Ac,P.Bc,P.C,[0,0],'InputDelay',[P.tau,0]);
P.sysd = c2d(P.sys,P.Ts);

%% Porownanie odpowiedzi modeli na skoki wejsc
outdir = fullfile(fileparts(mfilename('fullpath')),'wyniki_czesc1');
if ~isfolder(outdir), mkdir(outdir); end

% Dodatnie i ujemne skoki, od 10% do 50% doplywu w punkcie pracy.
steps = [-50 -25 -10 10 25 50];
n = 300; t = (0:n)'*P.Ts;
stepU = step(P.sysd(:,1),t); % odpowiedz na jednostkowy skok F1in
stepD = step(P.sysd(:,2),t); % odpowiedz na jednostkowy skok FD
fig=figure('Visible','off','Color','w'); hold on;
co=lines(numel(steps));
stepResults = zeros(numel(steps),5);
stepTracesU = zeros(numel(t),2*numel(steps));
for j=1:numel(steps)
    du=steps(j);
    u = (P.u0+du)*ones(n+1,1);
    d = P.d0*ones(n+1,1);
    ynl = simulate_nonlinear(P,u,d);
    ylin = P.h20 + du*stepU;
    hnl=((P.q0+du)/P.a2)^2;
    hlin=P.h20+du/P.k2;
    maxerr=max(abs(ynl-ylin));
    stepResults(j,:) = [du,hnl,hlin,hnl-hlin,maxerr];
    stepTracesU(:,2*j-1:2*j) = [ynl,ylin];
    plot(t,ynl,'-','Color',co(j,:),'LineWidth',1.4,'DisplayName',sprintf('%+d nieliniowy',du));
    plot(t,ylin,'--','Color',co(j,:),'LineWidth',1.2,'DisplayName',sprintf('%+d liniowy',du));
end
xlabel('t [s]'); ylabel('h_2 [cm]'); title('Skoki F_{1in} z punktu pracy');
grid on; legend('Location','eastoutside');
exportgraphics(fig,fullfile(outdir,'skoki_F1in.png'),'Resolution',160); close(fig);

% Zaklocenie dochodzi bez opoznienia transportowego sterowania.
fig=figure('Visible','off','Color','w'); hold on;
stepTracesD = zeros(numel(t),2*numel(steps));
distResults = zeros(numel(steps),5);
for j=1:numel(steps)
    dd=steps(j)*P.d0/P.u0; % te same procentowe zmiany modelu wzgłedem punktu pracy co w F1 - od 10% do 50% doplywu w punkcie pracy
    u = P.u0*ones(n+1,1);
    d = (P.d0+dd)*ones(n+1,1);
    ynl = simulate_nonlinear(P,u,d);
    ylin = P.h20 + dd*stepD;
    stepTracesD(:,2*j-1:2*j) = [ynl,ylin];
    hnl = ((P.q0+dd)/P.a2)^2;
    hlin = P.h20+dd/P.k2;
    distResults(j,:) = [dd,hnl,hlin,hnl-hlin,max(abs(ynl-ylin))];
    plot(t,ynl,'-','Color',co(j,:),'LineWidth',1.4,'DisplayName',sprintf('%+g nieliniowy',dd));
    plot(t,ylin,'--','Color',co(j,:),'LineWidth',1.2,'DisplayName',sprintf('%+g liniowy',dd));
end
xlabel('t [s]'); ylabel('h_2 [cm]'); title('Skoki F_D z punktu pracy');
grid on; legend('Location','eastoutside');
exportgraphics(fig,fullfile(outdir,'skoki_FD.png'),'Resolution',160); close(fig);

%% Konwencjonalny analityczny DMC
% Delta U = (M'' M + lambda I)^(-1) M'' (Yzad - Y0).
% Pierwszy element wektora Delta U jest stosowany do obiektu.
D=250; N=50; Nu=10; lambda=0.05;
% step zwraca probki od t=0; s(i) odpowiada chwili i*Ts.
response = step(P.sysd(:,1),(0:D+N)'*P.Ts);
s = response(2:end);
M=zeros(N,Nu);
for i=1:N
    for j=1:Nu
        if i>=j, M(i,j)=s(i-j+1); end
    end
end
K=(M'*M+lambda*eye(Nu))\M';
% Macierz wplywu dawnych zmian sterowania na przyszle wyjscie.
freeMatrix = zeros(N,D-1);
for i = 1:N
    freeMatrix(i,:) = s(i+(1:D-1))-s(1:D-1);
end

%% DMC: male i duze skoki, dodatnie i ujemne
% Kolumny: zmiana h2,zad [cm], zmiana FD [cm^3/s].
scenarios = [ 2,  10; ...
              4,  20; ... % poprzedni test
              8,  30; ...
             -8, -10];
tc = (0:400)'*P.Ts;
referencePeriod = tc>=400 & tc<3800;
disturbancePeriod = tc>=3800;
dmcResults = zeros(size(scenarios,1),6);
dmcTraces = struct([]);
fig = figure('Visible','off','Color','w','Position',[100 100 1050 1000]);
tiledlayout(4,2);

for j = 1:size(scenarios,1)
    deltaRef = scenarios(j,1);
    deltaDist = scenarios(j,2);
    ref = P.h20 + deltaRef*(tc>=400);
    dist = P.d0 + deltaDist*(tc>=3800);

    [yNL,uNL] = simulate_dmc(P,ref,dist,K,freeMatrix,'nonlinear');
    [yLin,uLin] = simulate_dmc(P,ref,dist,K,freeMatrix,'linear');
    dmcTraces(j).yNL = yNL;
    dmcTraces(j).uNL = uNL;
    dmcTraces(j).yLin = yLin;
    dmcTraces(j).uLin = uLin;
    dmcTraces(j).ref = ref;
    dmcTraces(j).dist = dist;

    iaeRef = P.Ts*sum(abs(yNL(referencePeriod)-ref(referencePeriod)));
    iaeDist = P.Ts*sum(abs(yNL(disturbancePeriod)-ref(disturbancePeriod)));
    maxDifference = max(abs(yNL-yLin));
    dmcResults(j,:) = [deltaRef,deltaDist,iaeRef,iaeDist, ...
                       maxDifference,min(uNL)];

    outputTile = 4*floor((j-1)/2)+mod(j-1,2)+1;
    nexttile(outputTile);
    plot(tc,yNL,'b',tc,yLin,'r--',tc,ref,'k:','LineWidth',1.2);
    grid on; ylabel('y(k) [cm]');
    title(sprintf('Skok y_{zad}: %+g cm, F_D: %+g cm^3/s',deltaRef,deltaDist));
    if j==1
        legend('obiekt nieliniowy','obiekt liniowy','y_{zad}(k)', ...
               'Location','best');
    end

    nexttile(outputTile+2);
    stairs(tc,[uNL;uNL(end)],'b','LineWidth',1.2);
    hold on;
    stairs(tc,[uLin;uLin(end)],'r--','LineWidth',1.2);
    grid on; xlabel('t [s]'); ylabel('u(k) [cm^3/s]');
    plot_dmc_scenario(tc,yNL,yLin,uNL,uLin,ref,deltaRef,deltaDist,j,outdir);

    if j==2
        % Zachowaj dotychczasowy wykres dla pierwotnego testu.
        plot_dmc_original(tc,yNL,uNL,ref,dist,outdir);
    end
end
exportgraphics(fig,fullfile(outdir,'regulacja_DMC_porownanie.png'), ...
               'Resolution',160);
close(fig);

stepSummary = array2table(stepResults,'VariableNames', ...
    {'Skok_F1in','h2_nieliniowy','h2_liniowy','Blad_ustalony','Maks_blad'});
disp(stepSummary);
dmcSummary = array2table(dmcResults,'VariableNames', ...
    {'Skok_zadania_cm','Skok_FD','IAE_zadanie','IAE_zaklocenie', ...
     'Maks_roznica_modeli_cm','Min_F1in'});
disp(dmcSummary);

% Dane do sprawozdania; zapis nie zmienia modelu ani algorytmu regulacji.
distSummary = array2table(distResults,'VariableNames', ...
    {'Skok_FD','h2_nieliniowy','h2_liniowy','Blad_ustalony','Maks_blad'});
writetable(stepSummary,fullfile(outdir,'porownanie_F1in.csv'));
writetable(distSummary,fullfile(outdir,'porownanie_FD.csv'));
writetable(dmcSummary,fullfile(outdir,'wskazniki_DMC.csv'));
save(fullfile(outdir,'dane_sprawozdanie.mat'),'P','D','N','Nu','lambda', ...
    't','tc','steps','scenarios','stepTracesU','stepTracesD','dmcTraces', ...
    'stepResults','distResults','dmcResults','s','M','K','freeMatrix');

%% Funkcje pomocnicze
function plot_dmc_scenario(t,yNL,yLin,uNL,uLin,ref,deltaRef,deltaDist,j,outdir)
    % Wyjscie i sterowanie sa nierozlaczna para dla kazdej proby.
    fig = figure('Visible','off','Color','w','Position',[100 100 1000 650]);
    tiledlayout(2,1);
    ay = nexttile;
    plot(t,yNL,'b',t,yLin,'r--','LineWidth',1.2);
    hold on; stairs(t,ref,'k:','LineWidth',1.2);
    grid on; ylabel('y(k) [cm]');
    title(sprintf('Skok y_{zad}: %+g cm, F_D: %+g cm^3/s',deltaRef,deltaDist));
    legend('obiekt nieliniowy','obiekt liniowy','y_{zad}(k)','Location','best');
    au = nexttile;
    stairs(t,[uNL;uNL(end)],'b','LineWidth',1.2);
    hold on; stairs(t,[uLin;uLin(end)],'r--','LineWidth',1.2);
    grid on; xlabel('t [s]'); ylabel('u(k) [cm^3/s]');
    linkaxes([ay,au],'x');
    stem = sprintf('DMC_scenariusz_%d',j);
    exportgraphics(fig,fullfile(outdir,[stem '.png']),'Resolution',160);
    exportgraphics(fig,fullfile(outdir,[stem '.pdf']),'ContentType','vector');
    close(fig);
end

function y = simulate_nonlinear(P,u,d)
    % Na kazdym odcinku Ts strumienie sa stale. Stanami sa objetosci.
    n = numel(u)-1;
    y = zeros(n+1,1);
    y(1) = P.h20;
    v = [P.v10;P.v20];
    for k = 1:n
        if k>P.delay, delayedU = u(k-P.delay); else, delayedU = P.u0; end
        v = volume_step(P,v,delayedU,d(k));
        y(k+1) = (max(v(2),0)/P.C2)^(1/3);
    end
end

function [y,u] = simulate_dmc(P,ref,dist,K,freeMatrix,plant)
    n = numel(ref)-1;
    y = zeros(n+1,1);
    u = zeros(n,1);
    y(1) = P.h20;
    v = [P.v10;P.v20];
    x = zeros(2,1);
    pastDu = zeros(size(freeMatrix,2),1);

    for k = 1:n
        % Przewidywanie bez nowych zmian sterowania (swobodna odpowiedz).
        free = y(k) + freeMatrix*pastDu;

        % Z analitycznego rozwiazania stosujemy tylko pierwszy przyrost.
        du = K(1,:)*(ref(k)-free);
        if k==1, previousU = P.u0; else, previousU = u(k-1); end
        u(k) = previousU + du;
        pastDu = [du;pastDu(1:end-1)];

        if k>P.delay, delayedU = u(k-P.delay); else, delayedU = P.u0; end
        if strcmp(plant,'nonlinear')
            v = volume_step(P,v,delayedU,dist(k));
            y(k+1) = (max(v(2),0)/P.C2)^(1/3);
        else
            deltaInput = [delayedU-P.u0;dist(k)-P.d0];
            x = P.sysd.A*x + P.sysd.B*deltaInput;
            y(k+1) = P.h20 + P.sysd.C*x;
        end
    end
end

function plot_dmc_original(t,y,u,ref,dist,outdir)
    fig = figure('Visible','off','Color','w');
    tiledlayout(2,1);
    nexttile;
    plot(t,y,'b',t,ref,'k--','LineWidth',1.4);
    grid on; ylabel('y(k) [cm]');
    legend('obiekt nieliniowy','y_{zad}(k)','Location','best');
    title('DMC: skok y_{zad} +4 cm, skok F_D +20 cm^3/s');
    nexttile;
    stairs(t(1:end-1),u,'b','LineWidth',1.3);
    hold on; yyaxis right;
    stairs(t,dist,'r--','LineWidth',1.1);
    grid on; xlabel('t [s]'); ylabel('F_D [cm^3/s]');
    yyaxis left; ylabel('u(k) [cm^3/s]');
    ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
    exportgraphics(fig,fullfile(outdir,'regulacja_DMC.png'),'Resolution',160);
    close(fig);
end

function v = volume_step(P,v,u,d)
    % MATLAB dobiera krok calkowania wewnatrz pojedynczego okresu Ts.
    [~,trajectory] = ode45(@(~,state) tank_rhs(P,state,u,d),[0 P.Ts],v);
    v = trajectory(end,:)';
end

function dv = tank_rhs(P,v,u,d)
    h1 = max(v(1),0)/P.A1;
    h2 = (max(v(2),0)/P.C2)^(1/3);
    f2 = P.a1*sqrt(h1);
    f3 = P.a2*sqrt(h2);
    dv = [u+d-f2; f2-f3];
end
