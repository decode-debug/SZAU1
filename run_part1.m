%% SZAU, zadanie 1, czesc 1: zbiorniki i konwencjonalny DMC
% Uruchom z folderu projektu: run_part1
% Wymaga Control System Toolbox (ss, c2d, step).
clear; close all; clc;

%% Parametry i punkt pracy
P = tank_parameters();
P.h10 = P.h0(1);
P.h20 = P.h0(2);
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
plots = ResultPlots(outdir);

% Dodatnie i ujemne skoki, od 10% do 50% doplywu w punkcie pracy.
steps = [-50 -25 -10 10 25 50];
distSteps = steps*P.d0/P.u0;
n = 300; t = (0:n)'*P.Ts;
stepU = step(P.sysd(:,1),t); % odpowiedz na jednostkowy skok F1in
stepD = step(P.sysd(:,2),t); % odpowiedz na jednostkowy skok FD (bez opoznienia)

[stepResults,stepTracesU] = evaluate_steps(P,steps,stepU,n,'u');
[distResults,stepTracesD] = evaluate_steps(P,distSteps,stepD,n,'d');
plots.steps(t,stepTracesU,steps,'F_{1in}','skoki_F1in.pdf');
plots.steps(t,stepTracesD,distSteps,'F_D','skoki_FD.pdf');

%% Konwencjonalny analityczny DMC
% Delta U = (M' M + lambda I)^(-1) M' (Yzad - Y0); stosowany jest pierwszy element.
D = 200; N = 50; Nu = 10; lambda = 0.05;
response = step(P.sysd(:,1),(0:D+N)'*P.Ts);
s = response(2:end);
[M,freeMatrix] = dmc_matrices(s,N,Nu,D);
K = (M'*M + lambda*eye(Nu))\M';

%% DMC: male i duze skoki, dodatnie i ujemne
% Kolumny: zmiana h2,zad [cm], zmiana FD [cm^3/s].
scenarios = [ 2,  10; ...
              4,  20; ...
              8,  30; ...
             -8, -10];
refStart = 400;
distStart = 3800;
tc = (0:400)'*P.Ts;
referencePeriod = tc>=refStart & tc<distStart;
disturbancePeriod = tc>=distStart;
dmcResults = zeros(size(scenarios,1),6);
dmcTraces = struct([]);

for j = 1:size(scenarios,1)
    deltaRef = scenarios(j,1);
    deltaDist = scenarios(j,2);
    ref = P.h20 + deltaRef*(tc>=refStart);
    dist = P.d0 + deltaDist*(tc>=distStart);

    [yNL,uNL] = simulate_dmc(P,ref,dist,K,freeMatrix,'nonlinear');
    [yLin,uLin] = simulate_dmc(P,ref,dist,K,freeMatrix,'linear');
    dmcTraces(j).yNL = yNL;
    dmcTraces(j).uNL = uNL;
    dmcTraces(j).yLin = yLin;
    dmcTraces(j).uLin = uLin;
    dmcTraces(j).ref = ref;
    dmcTraces(j).dist = dist;

    mseRef = mean((yNL(referencePeriod)-ref(referencePeriod)).^2);
    mseDist = mean((yNL(disturbancePeriod)-ref(disturbancePeriod)).^2);
    dmcResults(j,:) = [deltaRef,deltaDist,mseRef,mseDist,max(abs(yNL-yLin)),min(uNL)];
end
plots.dmc(tc,dmcTraces,scenarios,refStart,distStart);

stepSummary = array2table(stepResults,'VariableNames', ...
    {'Skok_F1in','h2_nieliniowy','h2_liniowy','Blad_ustalony','Maks_blad'});
distSummary = array2table(distResults,'VariableNames', ...
    {'Skok_FD','h2_nieliniowy','h2_liniowy','Blad_ustalony','Maks_blad'});
dmcSummary = array2table(dmcResults,'VariableNames', ...
    {'Skok_zadania_cm','Skok_FD','MSE_zadanie','MSE_zaklocenie', ...
     'Maks_roznica_modeli_cm','Min_F1in'});
disp(stepSummary);
disp(dmcSummary);

% Zapis wynikow symulacji.
writetable(stepSummary,fullfile(outdir,'porownanie_F1in.csv'));
writetable(distSummary,fullfile(outdir,'porownanie_FD.csv'));
writetable(dmcSummary,fullfile(outdir,'wskazniki_DMC.csv'));
save(fullfile(outdir,'dane_sprawozdanie.mat'),'P','D','N','Nu','lambda', ...
    't','tc','steps','scenarios','stepTracesU','stepTracesD','dmcTraces', ...
    'stepResults','distResults','dmcResults','s','M','K','freeMatrix');

%% Funkcje pomocnicze
function [results,traces] = evaluate_steps(P,amplitudes,stepResp,n,channel)
    results = zeros(numel(amplitudes),5);
    traces = zeros(n+1,2*numel(amplitudes));
    for j = 1:numel(amplitudes)
        amp = amplitudes(j);
        u = (P.u0 + amp*(channel=='u'))*ones(n+1,1);
        d = (P.d0 + amp*(channel=='d'))*ones(n+1,1);
        ynl = simulate_nonlinear(P,u,d);
        ylin = P.h20 + amp*stepResp;
        hnl = ((P.q0+amp)/P.a2)^2;
        hlin = P.h20 + amp/P.k2;
        results(j,:) = [amp,hnl,hlin,hnl-hlin,max(abs(ynl-ylin))];
        traces(:,2*j-1:2*j) = [ynl,ylin];
    end
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
        free = y(k) + freeMatrix*pastDu;
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

function v = volume_step(P,v,u,d)
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
