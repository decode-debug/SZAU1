%% SZAU, zadanie 1, czesc 3: numeryczny rozmyty regulator FGPC-SL
% Najpierw uruchom run_part2, potem run_part3.
% Wymaga Control System Toolbox i Optimization Toolbox (quadprog, lsqlin).
clear; close all; clc;

dataFile = fullfile('wyniki_czesc2','dane_czesc2.mat');
if ~isfile(dataFile)
    error('Brak %s. Uruchom najpierw run_part2.',dataFile);
end
load(dataFile,'P','modelTuned','ctrl','controlCases','tControl', ...
    'refStart','distStart','traces');
outdir = fullfile(fileparts(mfilename('fullpath')),'wyniki_czesc3');
plots = ResultPlots(outdir);

%% Nastawy numerycznego regulatora FGPC-SL
N_sl = 50;
Nu_sl = 10;
lambdaSL = [2.0, 1.0, 0.5, 0.2, 0.0];
selectedVariant = 3;  % wybrany wariant: lambda = 0.5
options = optimoptions('quadprog','Algorithm','active-set','MaxIterations',1000,'Display','off');

% Dodatkowe proby walidacyjne o posrednich amplitudach i przeciwnych znakach.
validationCases = [3,-8; -3,12; 6,-10; -6,15];
allCases = [controlCases;validationCases];
nProject = size(controlCases,1);
nVariants = numel(lambdaSL);
taskMask = tControl>=refStart & tControl<distStart;
distMask = tControl>=distStart;

metricNames = {'MSE_zadanie','MSE_zaklocenie','MSE_laczny','Min_u','Max_u','Maks_zmiana_u', ...
    'Suma_zmian_u','Cofniecie_cm','Przeregulowanie_cm','Maks_blad_po_FD_cm'};
resultNames = [{'Scenariusz','Wariant','N','Nu','Lambda'},metricNames, ...
    {'Sredni_czas_decyzji_s','Maks_czas_decyzji_s'}];
results = zeros(size(allCases,1)*nVariants,numel(resultNames));
allSL = struct([]);

for j = 1:size(allCases,1)
    if j<=nProject
        ref = traces(j).ref;
        dist = traces(j).dist;
    else
        ref = P.h0(2)+allCases(j,1)*(tControl>=refStart);
        dist = P.d0+allCases(j,2)*(tControl>=distStart);
    end
    for v = 1:nVariants
        cfg = ctrl;
        cfg.N = N_sl;
        cfg.lambda = lambdaSL(v);
        [y,u,decisionTime] = simulate_fgpc(P,modelTuned,cfg,Nu_sl,ref,dist,options,false);
        allSL(j,v).y = y;
        allSL(j,v).u = u;
        allSL(j,v).decisionTime = decisionTime;
        metrics = case_metrics(P,y,u,ref,taskMask,distMask);
        row = (j-1)*nVariants+v;
        results(row,:) = [j,v,N_sl,Nu_sl,lambdaSL(v),metrics, ...
            mean(decisionTime),max(decisionTime)];
    end
end

allScores = array2table(results,'VariableNames',resultNames);
mainMask = allScores.Scenariusz<=nProject;
slSummary = allScores(mainMask,:);
slTraces = allSL(1:nProject,:);
validationTraces = allSL(nProject+1:end,:);
validationSummary = allScores(~mainMask,:);
validationSummary.Scenariusz = validationSummary.Scenariusz-nProject;
validationSummary.Skok_zad_cm = validationCases(validationSummary.Scenariusz,1);
validationSummary.Skok_FD = validationCases(validationSummary.Scenariusz,2);

%% Zbiorcze porownanie: PDC (cz. 2) oraz FGPC-SL dla badanych wartosci lambda
methodNames = ["PDC";"PDC strojony"; ...
    "FGPC-SL (\lambda=2)";"FGPC-SL (\lambda=1)"; ...
    "FGPC-SL (\lambda=0,5)";"FGPC-SL (\lambda=0,2)";"FGPC-SL (\lambda=0)"];
comparison = zeros(nProject*numel(methodNames),2+numel(metricNames));
for j = 1:nProject
    first = (j-1)*numel(methodNames);
    for method = 1:2
        metrics = case_metrics(P,traces(j).y(:,method+2),traces(j).u(:,method+2), ...
            traces(j).ref,taskMask,distMask);
        comparison(first+method,:) = [j,method,metrics];
    end
    rows = slSummary.Scenariusz==j;
    comparison(first+(3:numel(methodNames)),:) = [slSummary.Scenariusz(rows), ...
        slSummary.Wariant(rows)+2,slSummary{rows,metricNames}];
end
comparisonSummary = array2table(comparison,'VariableNames',[{'Scenariusz','Regulator'},metricNames]);
comparisonSummary.Nazwa = methodNames(comparisonSummary.Regulator);
totalMSE = accumarray(comparisonSummary.Regulator, ...
    comparisonSummary.MSE_laczny,[],@mean);
totalRef = accumarray(comparisonSummary.Regulator,comparisonSummary.MSE_zadanie,[],@mean);
totalDist = accumarray(comparisonSummary.Regulator,comparisonSummary.MSE_zaklocenie,[],@mean);
disp(table(methodNames,totalRef,totalDist,totalMSE, ...
    'VariableNames',{'Regulator','MSE_zadanie','MSE_zaklocenie','MSE_laczny'}));

%% Wykresy i zapis wynikow
plots.plot_comparison(tControl,traces,slTraces,controlCases,lambdaSL,selectedVariant,refStart,distStart);
plots.plot_tuning(tControl,traces,slTraces,controlCases,lambdaSL,refStart,distStart);
writetable(slSummary,fullfile(outdir,'wskazniki_SL.csv'));
writetable(comparisonSummary,fullfile(outdir,'porownanie_regulatorow.csv'));
writetable(validationSummary,fullfile(outdir,'walidacja_SL.csv'));
save(fullfile(outdir,'dane_czesc3.mat'),'P','modelTuned','ctrl','controlCases', ...
    'tControl','traces','slTraces','N_sl','Nu_sl','lambdaSL','selectedVariant', ...
    'refStart','distStart','slSummary','comparisonSummary','totalMSE', ...
    'validationCases','validationTraces','validationSummary');

%% Funkcje pomocnicze
function metrics = case_metrics(P,y,u,ref,taskMask,distMask)
    err = y-ref;
    direction = sign(ref(end)-ref(1));
    level = direction*y(taskMask);
    retreat = max(cummax(level)-level);
    overshoot = max(0,max(direction*err(taskMask)));
    du = diff([P.u0;u]);
    metrics = [mean(err(taskMask).^2),mean(err(distMask).^2),mean(err(taskMask|distMask).^2), ...
        min(u),max(u),max(abs(du)),sum(abs(du)),retreat,overshoot,max(abs(err(distMask)))];
end
