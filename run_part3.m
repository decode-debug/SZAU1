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

%% Nastawy numerycznego regulatora FGPC-SL (wspolne N = 50, Nu = 10)
lambdaGrid = [2.0, 1.0, 0.5, 0.2, 0.0];
selectedVariant = 3;  % wybrany wariant: lambda = 0.5

% Dodatkowe proby walidacyjne o posrednich amplitudach i przeciwnych znakach.
validationCases = [3,-8; -3,12; 6,-10; -6,15];
allCases = [controlCases;validationCases];
nProject = size(controlCases,1);
nVariants = numel(lambdaGrid);
taskMask = tControl>=refStart & tControl<distStart;
distMask = tControl>=distStart;

metricNames = {'MSE_zadanie','MSE_zaklocenie','MSE_laczny','Przeregulowanie_cm', ...
    'Maks_blad_po_FD_cm','Min_u','Max_u','Maks_zmiana_u','Suma_zmian_u'};
resultNames = [{'Scenariusz','Wariant','N','Nu','Lambda'},metricNames, ...
    {'Sredni_czas_decyzji_s','Maks_czas_decyzji_s'}];
results = zeros(size(allCases,1)*nVariants,numel(resultNames));
allSL = struct([]);

for j = 1:size(allCases,1)
    ref = P.h0(2) + allCases(j,1)*(tControl>=refStart);
    dist = P.d0 + allCases(j,2)*(tControl>=distStart);
    for v = 1:nVariants
        cfg = ctrl;
        cfg.lambda = lambdaGrid(v);
        [y,u,dt] = simulate_fgpc(P,modelTuned,cfg,ref,dist,false);
        allSL(j,v).y = y;
        allSL(j,v).u = u;
        allSL(j,v).decisionTime = dt;
        m = control_metrics(P,y,u,ref,taskMask,distMask);
        row = (j-1)*nVariants+v;
        results(row,:) = [j,v,ctrl.N,ctrl.Nu,lambdaGrid(v),m,mean(dt),max(dt)];
    end
end

allScores = array2table(results,'VariableNames',resultNames);
mainMask = allScores.Scenariusz<=nProject;
slSummary = allScores(mainMask,:);
slTraces = allSL(1:nProject,:);
validationTraces = allSL(nProject+1:end,:);
validationSummary = allScores(~mainMask,:);
validationSummary.Scenariusz = validationSummary.Scenariusz - nProject;
validationSummary.Skok_zad_cm = validationCases(validationSummary.Scenariusz,1);
validationSummary.Skok_FD = validationCases(validationSummary.Scenariusz,2);

%% Zbiorcze porownanie: PDC (cz. 2) oraz FGPC-SL dla badanych wartosci lambda
methodNames = ["PDC";"PDC strojony"; ...
    arrayfun(@(l) sprintf("FGPC-SL (lambda=%g)",l),lambdaGrid.')];
comparison = zeros(nProject*numel(methodNames),2+numel(metricNames));
for j = 1:nProject
    first = (j-1)*numel(methodNames);
    for method = 1:2
        m = control_metrics(P,traces(j).y(:,method+2),traces(j).u(:,method+2), ...
            traces(j).ref,taskMask,distMask);
        comparison(first+method,:) = [j,method,m];
    end
    rows = slSummary.Scenariusz==j;
    comparison(first+(3:numel(methodNames)),:) = [slSummary.Scenariusz(rows), ...
        slSummary.Wariant(rows)+2,slSummary{rows,metricNames}];
end
comparisonSummary = array2table(comparison,'VariableNames',[{'Scenariusz','Regulator'},metricNames]);
comparisonSummary.Nazwa = methodNames(comparisonSummary.Regulator);
totalRef = accumarray(comparisonSummary.Regulator,comparisonSummary.MSE_zadanie,[],@mean);
totalDist = accumarray(comparisonSummary.Regulator,comparisonSummary.MSE_zaklocenie,[],@mean);
totalMSE = accumarray(comparisonSummary.Regulator,comparisonSummary.MSE_laczny,[],@mean);
disp(table(methodNames,totalRef,totalDist,totalMSE, ...
    'VariableNames',{'Regulator','MSE_zadanie','MSE_zaklocenie','MSE_laczny'}));

%% Wykresy i zapis wynikow
plots.plot_comparison(tControl,traces,slTraces,controlCases,lambdaGrid,selectedVariant,refStart,distStart);
plots.plot_tuning(tControl,traces,slTraces,controlCases,lambdaGrid,refStart,distStart,'FGPC-SL','strojenie_SL.pdf');
writetable(slSummary,fullfile(outdir,'wskazniki_SL.csv'));
writetable(comparisonSummary,fullfile(outdir,'porownanie_regulatorow.csv'));
writetable(validationSummary,fullfile(outdir,'walidacja_SL.csv'));
save(fullfile(outdir,'dane_czesc3.mat'),'P','modelTuned','ctrl','controlCases', ...
    'tControl','traces','slTraces','lambdaGrid','selectedVariant', ...
    'refStart','distStart','slSummary','comparisonSummary','totalMSE', ...
    'validationCases','validationTraces','validationSummary');
