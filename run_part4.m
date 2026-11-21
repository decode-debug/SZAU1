%% SZAU, zadanie 1, czesc 4: numeryczny rozmyty regulator FGPC-NPL
% Najpierw uruchom run_part2 i run_part3, potem run_part4.
% Algorytm zgodny ze skryptem SZAU (s. 168-169, wzory 4.103-4.105 i 4.118-4.125).
clear; close all; clc;

root = fileparts(mfilename('fullpath'));
source = fullfile(root,'wyniki_czesc3','dane_czesc3.mat');
if ~isfile(source), error('Uruchom najpierw run_part2 i run_part3.'); end
load(source,'P','ctrl','modelTuned','traces','slTraces','controlCases', ...
    'tControl','refStart','distStart','N_sl','Nu_sl','lambdaSL','selectedVariant', ...
    'validationCases','validationTraces');
outdir = fullfile(root,'wyniki_czesc4');
plots = ResultPlots(outdir);
options = optimoptions('quadprog','Algorithm','active-set','Display','off');

%% Symulacje FGPC-NPL dla tych samych nastaw N, Nu i zestawu lambda co FGPC-SL
% Porownanie przy tej samej lambda oraz wplyw strojenia lambda (Przyklad 4.19, Rys. 4.44-4.45).
N_npl = N_sl;
Nu_npl = Nu_sl;
lambdaNPL = lambdaSL;  % [2.0, 1.0, 0.5, 0.2]
nProject = numel(traces);
nVariants = numel(lambdaNPL);
allCases = [controlCases;validationCases];

taskMask = tControl>=refStart & tControl<distStart;
distMask = tControl>=distStart;

allNPL = struct([]);
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
        cfg.N = N_npl;
        cfg.lambda = lambdaNPL(v);
        [y,u,dt] = simulate_fgpc(P,modelTuned,cfg,Nu_npl,ref,dist,options,true);
        allNPL(j,v).y = y;
        allNPL(j,v).u = u;
        allNPL(j,v).decisionTime = dt;
    end
end
nplTraces = allNPL(1:nProject,:);
nplValTraces = allNPL(nProject+1:end,:);

%% Zbiorcze zestawienie wszystkich regulatorow z czesci 1-4
% Wybrane warianty: DMC cz.1, DMC ogr., PDC, PDC strojony,
% FGPC-SL (lam=1, lam=0.5, lam=0.2) oraz FGPC-NPL (lam=1, lam=0.5, lam=0.2).
names = [ ...
    "DMC cz. 1"; ...
    "DMC ograniczony"; ...
    "PDC"; ...
    "PDC strojony"; ...
    "FGPC-SL (lambda=1)"; ...
    "FGPC-SL (lambda=0,5)"; ...
    "FGPC-SL (lambda=0,2)"; ...
    "FGPC-NPL (lambda=1)"; ...
    "FGPC-NPL (lambda=0,5)"; ...
    "FGPC-NPL (lambda=0,2)"];

nMethods = numel(names);
rows = zeros(nProject*nMethods,11);
for j = 1:nProject
    ref = traces(j).ref;
    Y = [traces(j).y, ...
         slTraces(j,2).y, slTraces(j,3).y, slTraces(j,4).y, ...
         nplTraces(j,2).y, nplTraces(j,3).y, nplTraces(j,4).y];
    U = [traces(j).u, ...
         slTraces(j,2).u, slTraces(j,3).u, slTraces(j,4).u, ...
         nplTraces(j,2).u, nplTraces(j,3).u, nplTraces(j,4).u];
    for m = 1:nMethods
        e = Y(:,m)-ref;
        du = diff([P.u0;U(:,m)]);
        dir = sign(controlCases(j,1));
        z = dir*Y(taskMask,m);
        over = max(0,max(dir*e(taskMask)));
        retr = max(cummax(z)-z);
        rows((j-1)*nMethods+m,:) = [j,m, ...
            P.Ts*sum(abs(e(taskMask))), P.Ts*sum(abs(e(distMask))), ...
            over, retr, max(abs(e(distMask))), ...
            min(U(:,m)), max(U(:,m)), max(abs(du)), sum(abs(du))];
    end
end

quality = array2table(rows,'VariableNames',{ ...
    'Scenariusz','ID','IAE_zadanie','IAE_zaklocenie', ...
    'Przeregulowanie_cm','Cofniecie_cm','Maks_blad_po_FD_cm', ...
    'Min_u','Max_u','Maks_zmiana_u','Suma_zmian_u'});
quality.Regulator = names(quality.ID);

totals = table(names, ...
    accumarray(quality.ID,quality.IAE_zadanie), ...
    accumarray(quality.ID,quality.IAE_zaklocenie), ...
    accumarray(quality.ID,quality.Przeregulowanie_cm,[],@max), ...
    accumarray(quality.ID,quality.Cofniecie_cm,[],@max), ...
    accumarray(quality.ID,quality.Maks_blad_po_FD_cm,[],@max), ...
    accumarray(quality.ID,quality.Suma_zmian_u), ...
    'VariableNames',{'Regulator','IAE_zadanie','IAE_zaklocenie', ...
    'Maks_przereg_cm','Maks_cofniecie_cm','Maks_blad_FD_cm','Suma_zmian_u'});
totals.IAE = totals.IAE_zadanie + totals.IAE_zaklocenie;
disp(totals);

% Zestawienie wplywu lambda dla SL i NPL (wszystkie 4 wartosci lambdaSL).
sweepRows = zeros(nVariants,11);
for v = 1:nVariants
    slR=0; slD=0; slO=0; slRt=0; slDU=0;
    npR=0; npD=0; npO=0; npRt=0; npDU=0;
    for j = 1:nProject
        ref = traces(j).ref; dir = sign(controlCases(j,1));
        eS = slTraces(j,v).y-ref; zS = dir*slTraces(j,v).y(taskMask);
        slR=slR+P.Ts*sum(abs(eS(taskMask))); slD=slD+P.Ts*sum(abs(eS(distMask)));
        slO=max(slO,max(0,max(dir*eS(taskMask)))); slRt=max(slRt,max(cummax(zS)-zS));
        slDU=slDU+sum(abs(diff([P.u0;slTraces(j,v).u])));

        eN = nplTraces(j,v).y-ref; zN = dir*nplTraces(j,v).y(taskMask);
        npR=npR+P.Ts*sum(abs(eN(taskMask))); npD=npD+P.Ts*sum(abs(eN(distMask)));
        npO=max(npO,max(0,max(dir*eN(taskMask)))); npRt=max(npRt,max(cummax(zN)-zN));
        npDU=npDU+sum(abs(diff([P.u0;nplTraces(j,v).u])));
    end
    sweepRows(v,:) = [lambdaNPL(v), slR, slD, slR+slD, slO, slDU, ...
                                    npR, npD, npR+npD, npO, npDU];
end
sweepTable = array2table(sweepRows,'VariableNames',{ ...
    'Lambda','SL_IAE_zad','SL_IAE_zakl','SL_IAE','SL_Przereg_cm','SL_Suma_du', ...
    'NPL_IAE_zad','NPL_IAE_zakl','NPL_IAE','NPL_Przereg_cm','NPL_Suma_du'});
disp(sweepTable);

% Walidacja NPL vs SL na 4 probach dodatkowych (dla wybranej lam=0.5).
valRows = zeros(size(validationCases,1),8);
for j = 1:size(validationCases,1)
    ref = P.h0(2)+validationCases(j,1)*(tControl>=refStart);
    eS = validationTraces(j,selectedVariant).y-ref;
    eN = nplValTraces(j,selectedVariant).y-ref;
    valRows(j,:) = [validationCases(j,1), validationCases(j,2), ...
        P.Ts*sum(abs(eS(taskMask))), P.Ts*sum(abs(eS(distMask))), P.Ts*sum(abs(eS(taskMask|distMask))), ...
        P.Ts*sum(abs(eN(taskMask))), P.Ts*sum(abs(eN(distMask))), P.Ts*sum(abs(eN(taskMask|distMask)))];
end
valTable = array2table(valRows,'VariableNames',{ ...
    'Skok_zad_cm','Skok_FD','SL_IAE_zad','SL_IAE_zakl','SL_IAE', ...
    'NPL_IAE_zad','NPL_IAE_zakl','NPL_IAE'});
disp(valTable);

writetable(quality,fullfile(outdir,'porownanie_NPL.csv'));
writetable(totals,fullfile(outdir,'podsumowanie_NPL.csv'));
writetable(sweepTable,fullfile(outdir,'strojenie_NPL.csv'));
writetable(valTable,fullfile(outdir,'walidacja_NPL.csv'));

results = nplTraces(:,selectedVariant);  % glowny wariant NPL (lambda = 0.5)
save(fullfile(outdir,'dane_czesc4.mat'),'results','nplTraces','nplValTraces', ...
    'tControl','controlCases','validationCases','ctrl','N_npl','Nu_npl', ...
    'lambdaNPL','selectedVariant','quality','totals','sweepTable','valTable');

plots.npl(tControl,traces,slTraces,nplTraces,selectedVariant,lambdaNPL,controlCases,refStart,distStart);
