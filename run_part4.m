%% SZAU, zadanie 1, czesc 4: numeryczny rozmyty regulator FGPC-NPL
% Najpierw uruchom run_part2 i run_part3, potem run_part4.
clear; close all; clc;

root = fileparts(mfilename('fullpath'));
source = fullfile(root,'wyniki_czesc3','dane_czesc3.mat');
if ~isfile(source), error('Uruchom najpierw run_part2 i run_part3.'); end
load(source,'P','ctrl','modelTuned','traces','slTraces','controlCases', ...
    'tControl','refStart','distStart','lambdaGrid','selectedVariant', ...
    'validationCases','validationTraces');
outdir = fullfile(root,'wyniki_czesc4');
plots = ResultPlots(outdir);

%% Symulacje FGPC-NPL dla tych samych nastaw N, Nu i zestawu lambda co FGPC-SL
nProject = numel(traces);
nVariants = numel(lambdaGrid);
allCases = [controlCases;validationCases];
taskMask = tControl>=refStart & tControl<distStart;
distMask = tControl>=distStart;

allNPL = struct([]);
for j = 1:size(allCases,1)
    ref = P.h0(2) + allCases(j,1)*(tControl>=refStart);
    dist = P.d0 + allCases(j,2)*(tControl>=distStart);
    for v = 1:nVariants
        cfg = ctrl;
        cfg.lambda = lambdaGrid(v);
        [y,u,dt] = simulate_fgpc(P,modelTuned,cfg,ref,dist,true);
        allNPL(j,v).y = y;
        allNPL(j,v).u = u;
        allNPL(j,v).decisionTime = dt;
    end
end
nplTraces = allNPL(1:nProject,:);
nplValTraces = allNPL(nProject+1:end,:);

%% Zbiorcze zestawienie regulatorow z czesci 1-4
varIdx = 2:nVariants; % lambda = 1, 0.5, 0.2, 0
names = ["DMC cz. 1"; "DMC ograniczony"; "PDC"; "PDC strojony"; ...
    arrayfun(@(l) sprintf("FGPC-SL (lambda=%g)",l),lambdaGrid(varIdx).'); ...
    arrayfun(@(l) sprintf("FGPC-NPL (lambda=%g)",l),lambdaGrid(varIdx).')];

nMethods = numel(names);
rows = zeros(nProject*nMethods,11);
for j = 1:nProject
    ref = traces(j).ref;
    Y = [traces(j).y, slTraces(j,varIdx).y, nplTraces(j,varIdx).y];
    U = [traces(j).u, slTraces(j,varIdx).u, nplTraces(j,varIdx).u];
    for m = 1:nMethods
        rows((j-1)*nMethods+m,:) = [j, m, control_metrics(P,Y(:,m),U(:,m),ref,taskMask,distMask)];
    end
end

quality = array2table(rows,'VariableNames',{ ...
    'Scenariusz','ID','MSE_zadanie','MSE_zaklocenie','MSE_laczny', ...
    'Przeregulowanie_cm','Maks_blad_po_FD_cm','Min_u','Max_u','Maks_zmiana_u','Suma_zmian_u'});
quality.Regulator = names(quality.ID);

totals = table(names, ...
    accumarray(quality.ID,quality.MSE_zadanie,[],@mean), ...
    accumarray(quality.ID,quality.MSE_zaklocenie,[],@mean), ...
    accumarray(quality.ID,quality.MSE_laczny,[],@mean), ...
    accumarray(quality.ID,quality.Przeregulowanie_cm,[],@max), ...
    accumarray(quality.ID,quality.Maks_blad_po_FD_cm,[],@max), ...
    accumarray(quality.ID,quality.Suma_zmian_u), ...
    'VariableNames',{'Regulator','MSE_zadanie','MSE_zaklocenie','MSE_laczny', ...
    'Maks_przereg_cm','Maks_blad_FD_cm','Suma_zmian_u'});
disp(totals);

% Zestawienie wplywu lambda dla SL i NPL.
sweepRows = zeros(nVariants,11);
for v = 1:nVariants
    mSL = zeros(nProject,9);
    mNPL = zeros(nProject,9);
    for j = 1:nProject
        ref = traces(j).ref;
        mSL(j,:)  = control_metrics(P,slTraces(j,v).y, slTraces(j,v).u, ref,taskMask,distMask);
        mNPL(j,:) = control_metrics(P,nplTraces(j,v).y,nplTraces(j,v).u,ref,taskMask,distMask);
    end
    sweepRows(v,:) = [lambdaGrid(v), ...
        mean(mSL(:,1:3),1),  max(mSL(:,4)),  sum(mSL(:,9)), ...
        mean(mNPL(:,1:3),1), max(mNPL(:,4)), sum(mNPL(:,9))];
end
sweepTable = array2table(sweepRows,'VariableNames',{ ...
    'Lambda','SL_MSE_zad','SL_MSE_zakl','SL_MSE','SL_Przereg_cm','SL_Suma_du', ...
    'NPL_MSE_zad','NPL_MSE_zakl','NPL_MSE','NPL_Przereg_cm','NPL_Suma_du'});
disp(sweepTable);

% Walidacja NPL vs SL na 4 probach dodatkowych (dla wybranej lam=0.5).
valRows = zeros(size(validationCases,1),8);
for j = 1:size(validationCases,1)
    ref = P.h0(2) + validationCases(j,1)*(tControl>=refStart);
    mS = control_metrics(P,validationTraces(j,selectedVariant).y, ...
        validationTraces(j,selectedVariant).u,ref,taskMask,distMask);
    mN = control_metrics(P,nplValTraces(j,selectedVariant).y, ...
        nplValTraces(j,selectedVariant).u,ref,taskMask,distMask);
    valRows(j,:) = [validationCases(j,1), validationCases(j,2), mS(1:3), mN(1:3)];
end
valTable = array2table(valRows,'VariableNames',{ ...
    'Skok_zad_cm','Skok_FD','SL_MSE_zad','SL_MSE_zakl','SL_MSE', ...
    'NPL_MSE_zad','NPL_MSE_zakl','NPL_MSE'});
disp(valTable);

writetable(quality,fullfile(outdir,'porownanie_NPL.csv'));
writetable(totals,fullfile(outdir,'podsumowanie_NPL.csv'));
writetable(sweepTable,fullfile(outdir,'strojenie_NPL.csv'));
writetable(valTable,fullfile(outdir,'walidacja_NPL.csv'));

results = nplTraces(:,selectedVariant);  % glowny wariant NPL (lambda = 0.5)
save(fullfile(outdir,'dane_czesc4.mat'),'results','nplTraces','nplValTraces', ...
    'tControl','controlCases','validationCases','ctrl', ...
    'lambdaGrid','selectedVariant','quality','totals','sweepTable','valTable');

plots.npl(tControl,traces,slTraces,nplTraces,selectedVariant,lambdaGrid,controlCases,refStart,distStart);
