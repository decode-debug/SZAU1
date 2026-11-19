%% SZAU, zadanie 1, czesc 3: numeryczny rozmyty DMC-SL
% Najpierw uruchom run_part2, potem run_part3.
% Wymaga Control System Toolbox, Fuzzy Logic Toolbox i Optimization Toolbox.
clear; close all; clc;

dataFile = fullfile('wyniki_czesc2','dane_czesc2.mat');
if ~isfile(dataFile)
    error('Brak %s. Uruchom najpierw run_part2.',dataFile);
end
load(dataFile,'P','modelTuned','ctrl','controlCases','tControl','traces');
outdir = fullfile(fileparts(mfilename('fullpath')),'wyniki_czesc3');
plots = ResultPlots(outdir);

% Kolejno SL o nastawach PDC, poprzedni wariant i wariant po strojeniu.
predictionHorizons = [ctrl.N,50,20];
horizons = [1,5,3];
lambdaSL = [ctrl.lambda,0.2,0.005];
referenceTime = [0,0,160];   % [s] stala czasowa filtru wartosci zadanej
localSteps = local_step_responses(P,modelTuned,ctrl.D);
options = optimoptions('quadprog','Algorithm','active-set','Display','off');

% Dodatkowe amplitudy i przeciwne znaki zadania i zaklocenia.
validationCases = [3,-8; -3,12; 6,-10; -6,15];
allCases = [controlCases;validationCases];
nProject = size(controlCases,1);
nVariants = numel(horizons);
refStart = 400; distStart = 3200;
taskMask = tControl>=refStart & tControl<distStart;
distMask = tControl>=distStart;
metricNames = {'IAE_zadanie','IAE_zaklocenie','Min_u','Max_u','Maks_zmiana_u', ...
    'Suma_zmian_u','Cofniecie_cm','Przeregulowanie_cm','Maks_blad_po_FD_cm'};
resultNames = [{'Scenariusz','Wariant','N','Nu','Lambda','Tf_s'},metricNames, ...
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
        cfg.N = predictionHorizons(v);
        cfg.lambda = lambdaSL(v);
        [y,u,decisionTime] = simulate_fdmc(P,localSteps,modelTuned.centers, ...
            cfg,horizons(v),referenceTime(v),ref,dist,options);
        allSL(j,v).y = y;
        allSL(j,v).u = u;
        allSL(j,v).decisionTime = decisionTime;
        metrics = case_metrics(P,y,u,ref,taskMask,distMask);
        row = (j-1)*nVariants+v;
        results(row,:) = [j,v,cfg.N,horizons(v),cfg.lambda,referenceTime(v), ...
            metrics,mean(decisionTime),max(decisionTime)];
    end
end
allScores = array2table(results,'VariableNames',resultNames);
main = allScores.Scenariusz<=nProject;
slSummary = allScores(main,:);
slTraces = allSL(1:nProject,:);
validationTraces = allSL(nProject+1:end,:);
validationSummary = allScores(~main,:);
validationSummary.Scenariusz = validationSummary.Scenariusz-nProject;
validationSummary.Skok_zad_cm = validationCases(validationSummary.Scenariusz,1);
validationSummary.Skok_FD = validationCases(validationSummary.Scenariusz,2);

% Identyfikatory 1 i 2 to analityczne PDC z czesci 2, kolejne dotycza SL.
methodNames = ["PDC";"PDC strojony";"SL o nastawach PDC";"SL przed strojeniem";"SL po strojeniu"];
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
totalIAE = accumarray(comparisonSummary.Regulator, ...
    comparisonSummary.IAE_zadanie+comparisonSummary.IAE_zaklocenie);
disp(table(methodNames,totalIAE,'VariableNames',{'Regulator','Suma_IAE'}));

plots.plot_comparison(tControl,traces,slTraces,controlCases,methodNames,refStart,distStart);
plots.plot_tuning(tControl,traces,slTraces,controlCases,refStart,distStart);
writetable(slSummary,fullfile(outdir,'wskazniki_SL.csv'));
writetable(comparisonSummary,fullfile(outdir,'porownanie_regulatorow.csv'));
writetable(validationSummary,fullfile(outdir,'walidacja_SL.csv'));
save(fullfile(outdir,'dane_czesc3.mat'),'P','modelTuned','ctrl','controlCases', ...
    'tControl','traces','slTraces','horizons','predictionHorizons','lambdaSL', ...
    'referenceTime','refStart','distStart','localSteps','slSummary','comparisonSummary', ...
    'totalIAE','validationCases','validationTraces','validationSummary');

%% Funkcje pomocnicze
function metrics = case_metrics(P,y,u,ref,taskMask,distMask)
    err = y-ref;
    direction = sign(ref(end)-ref(1));
    level = direction*y(taskMask);
    retreat = max(cummax(level)-level);
    overshoot = max(0,max(direction*err(taskMask)));
    du = diff([P.u0;u]);
    metrics = [P.Ts*sum(abs(err(taskMask))),P.Ts*sum(abs(err(distMask))), ...
        min(u),max(u),max(abs(du)),sum(abs(du)),retreat,overshoot,max(abs(err(distMask)))];
end
