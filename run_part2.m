%% SZAU, zadanie 1, czesc 2: dyskretne modele Takagiego-Sugeno i rozmyty analityczny DMC (PDC)
% Uruchom z folderu projektu: run_part2
% Wymaga Control System Toolbox (ss, c2d, step, tfdata).
clear; close all; clc;

%% Parametry obiektu (identyczne jak w czesci 1)
P = tank_parameters();

outdir = fullfile(fileparts(mfilename('fullpath')),'wyniki_czesc2');
plots = ResultPlots(outdir);

%% Modele rozmyte TS (wzory 4.38 i 4.41): od 2 do 5 modeli lokalnych
% Zmienna przeslanek: mierzony poziom h2 (silnie zmienia pojemnosc 3*C2*h2^2
% i odzwierciedla wplyw niemierzalnego zaklocenia FD).
% Zakres odpowiada zmianom F1in = 50..150 cm^3/s (+/-50%) przy FD = 20 cm^3/s.
hRange = [(70/P.a2)^2, (170/P.a2)^2];
modelCounts = 2:5;
models = cell(size(modelCounts));
for j = 1:numel(modelCounts)
    centers = linspace(hRange(1),hRange(2),modelCounts(j));
    models{j} = build_arx_model(P,centers);
end
linearModel = build_arx_model(P,P.h0(2));

%% Ocena modeli: 12 skokow F1in i FD (+/-10%, +/-25%, +/-50%) z punktu pracy
openSteps = [-50 0; -25 0; -10 0; 10 0; 25 0; 50 0; ...
              0 -10; 0 -5; 0 -2; 0 2; 0 5; 0 10];
nOpen = 300;
tOpen = (0:nOpen)'*P.Ts;
truth = zeros(nOpen+1,size(openSteps,1));
for j = 1:size(openSteps,1)
    truth(:,j) = simulate_plant_step(P,openSteps(j,:),nOpen);
end

modelRmse = zeros(numel(models),1);
modelMaxError = zeros(numel(models),1);
stepMse = zeros(size(openSteps,1),numel(models)+2);
modelOutputs = zeros([size(truth),numel(models)]);
linearOutputs = zeros(size(truth));
linearErrors = zeros(size(truth));
for j = 1:size(openSteps,1)
    linearOutputs(:,j) = simulate_arx_step(P,linearModel,openSteps(j,:),nOpen);
    linearErrors(:,j) = linearOutputs(:,j)-truth(:,j);
end
linearRmse = sqrt(mean(linearErrors(:).^2));
linearMaxError = max(abs(linearErrors(:)));
stepMse(:,1) = mean(linearErrors.^2,1).';

for m = 1:numel(models)
    errors = zeros(size(truth));
    for j = 1:size(openSteps,1)
        yTS = simulate_arx_step(P,models{m},openSteps(j,:),nOpen);
        modelOutputs(:,j,m) = yTS;
        errors(:,j) = yTS-truth(:,j);
    end
    modelRmse(m) = sqrt(mean(errors(:).^2));
    modelMaxError(m) = max(abs(errors(:)));
    stepMse(:,m+1) = mean(errors.^2,1).';
end
modelSummary = table([1;modelCounts'],[linearRmse;modelRmse], ...
    [linearMaxError;modelMaxError], ...
    'VariableNames',{'Liczba_regul','RMSE_cm','Maks_blad_cm'});
disp(modelSummary);

%% Wybor modelu (5 regul - najmniejszy blad) i strojenie funkcji przynaleznosci
selectedCount = 5;
selectedIndex = find(modelCounts==selectedCount);
modelInitial = models{selectedIndex};

% Przesuwamy trzy wewnetrzne wezly funkcji przynaleznosci, pozostawiajac
% konce zakresu i lokalne nastepniki stale.
middle = modelInitial.centers(2:4);
shifts = [-2 -1 0 1 2];
bestRmse = inf;
modelTuned = modelInitial;
tunedOutputs = modelOutputs(:,:,selectedIndex);
candidateOutputs = zeros(size(truth));
for s1 = shifts
    for s2 = shifts
        for s3 = shifts
            candidateCenters = [hRange(1),middle+[s1,s2,s3],hRange(2)];
            if any(diff(candidateCenters)<=2), continue; end
            candidateModel = modelInitial;
            candidateModel.centers = candidateCenters;
            candidateMse = zeros(size(openSteps,1),1);
            for j = 1:size(openSteps,1)
                yTS = simulate_arx_step(P,candidateModel,openSteps(j,:),nOpen);
                err = yTS-truth(:,j);
                candidateOutputs(:,j) = yTS;
                candidateMse(j) = mean(err.^2);
            end
            rmse = sqrt(mean(candidateMse));
            if rmse<bestRmse
                bestRmse = rmse;
                modelTuned = candidateModel;
                tunedOutputs = candidateOutputs;
                stepMse(:,end) = candidateMse;
            end
        end
    end
end
fprintf('Model TS 5: RMSE przed = %.4f cm, po = %.4f cm; wezly = %s cm\n', ...
    modelRmse(selectedIndex),bestRmse,mat2str(modelTuned.centers,4));

% Walidacja modelu na 4 dodatkowych skokach nieuzytych przy strojeniu.
validationSteps = [-40 0; 40 0; 0 -7; 0 7];
validationError = zeros(nOpen+1,size(validationSteps,1),2);
for j = 1:size(validationSteps,1)
    yTrue = simulate_plant_step(P,validationSteps(j,:),nOpen);
    validationError(:,j,1) = simulate_arx_step(P,modelInitial,validationSteps(j,:),nOpen)-yTrue;
    validationError(:,j,2) = simulate_arx_step(P,modelTuned,validationSteps(j,:),nOpen)-yTrue;
end
validationRmse = [sqrt(mean(validationError(:,:,1).^2,'all')); ...
                  sqrt(mean(validationError(:,:,2).^2,'all'))];
fprintf('Walidacja TS 5: RMSE przed = %.4f cm, po = %.4f cm\n', ...
    validationRmse(1),validationRmse(2));

% Podsumowanie wg amplitudy skoku (10%, 25%, 50%).
amplitudePct = 100*max(abs(openSteps)./[P.u0,P.d0],[],2);
ampLevels = unique(amplitudePct);
groupRmse = zeros(numel(ampLevels),size(stepMse,2));
for j = 1:numel(ampLevels)
    groupRmse(j,:) = sqrt(mean(stepMse(amplitudePct==ampLevels(j),:),1));
end
amplitudeSummary = array2table([ampLevels,groupRmse],'VariableNames', ...
    {'Skok_proc','Liniowy_cm','TS2_cm','TS3_cm','TS4_cm','TS5_cm','TS5_strojony_cm'});
disp(amplitudeSummary);

% Dryf modeli przy zerowym wymuszeniu po jednym czasie opoznienia.
driftModels = [{linearModel},models,{modelTuned}];
nominalDrift = zeros(1,numel(driftModels));
for j = 1:numel(driftModels)
    yNom = simulate_arx_step(P,driftModels{j},[0,0],P.delay);
    nominalDrift(j) = yNom(end)-P.h0(2);
end

%% Wykresy modeli i funkcji przynaleznosci
levels = linspace(hRange(1),hRange(2),500)';
modelWeights = cellfun(@(m) ts_memberships(levels,m.centers),models,'UniformOutput',false);
initialWeights = ts_memberships(levels,modelTuned.centers);
plots.plot_memberships(levels,modelWeights,modelCounts,initialWeights);
plots.plot_model_responses(modelOutputs(:,:,selectedIndex),tunedOutputs,openSteps,truth,tOpen);
plots.plot_all_models(linearOutputs,modelOutputs,modelCounts,openSteps,truth,tOpen);

%% Analityczny rozmyty DMC (PDC, wzory 4.72-4.75) z ograniczeniami sterowania
% Te same nastawy bazowe co konwencjonalny DMC z czesci 1: D=250, N=50, Nu=10, lambda=0.05.
% Ograniczenia sterowania: 0 <= u <= 200 cm^3/s oraz |du| <= 15 cm^3/s na probke.
ctrl.D = 250;
ctrl.N = 50;
ctrl.Nu = 10;
ctrl.lambda = 0.05;
ctrl.uMin = 0;         % cm^3/s
ctrl.uMax = 200;       % cm^3/s
ctrl.duMax = 15;       % cm^3/s na probke

linearControllers = design_dmc(P,linearModel,ctrl);
fuzzyControllers = design_dmc(P,modelTuned,ctrl);

% Ta sama os czasu i te same 4 scenariusze co w czesci 1 (skok zadanej 400 s, skok FD 3800 s).
controlCases = [ 2,  10; ...
                 4,  20; ...
                 8,  30; ...
                -8, -10];
refStart = 400;
distStart = 3800;
tControl = (0:400)'*P.Ts;
refMask = tControl>=refStart & tControl<distStart;
distMask = tControl>=distStart;

testData = cell(size(controlCases,1),1);
for j = 1:numel(testData)
    testData{j}.ref = P.h0(2)+controlCases(j,1)*(tControl>=refStart);
    testData{j}.dist = P.d0+controlCases(j,2)*(tControl>=distStart);
end

% Strojenie funkcji przynaleznosci regulatora PDC (wykladnik gamma + przesuniecia wezlow).
gammaInitial = 1;
gammaCandidates = [0.5 0.75 1 1.25 1.5 2 3];
gammaScores = zeros(size(gammaCandidates));
for g = 1:numel(gammaCandidates)
    gammaScores(g) = pdc_mse(P,fuzzyControllers,modelTuned.centers, ...
        gammaCandidates(g),testData,ctrl,refMask,distMask);
end
[~,bestGammaIndex] = min(gammaScores);
gammaTuned = gammaCandidates(bestGammaIndex);
controllerCenters = modelTuned.centers;
bestControlScore = gammaScores(bestGammaIndex);

% Przesuniecia wewnetrznych wezlow w regulatorze (stale lokalne regulatory).
cShifts = [-2 0 2];
for d1 = cShifts
    for d2 = cShifts
        for d3 = cShifts
            candCenters = modelTuned.centers;
            candCenters(2:4) = candCenters(2:4) + [d1,d2,d3];
            if any(diff(candCenters)<=1.5), continue; end
            score = pdc_mse(P,fuzzyControllers,candCenters,gammaTuned,testData,ctrl,refMask,distMask);
            if score<bestControlScore
                bestControlScore = score;
                controllerCenters = candCenters;
            end
        end
    end
end

%% Symulacje porownawcze 4 regulatorow czesci 2
controlResults = zeros(size(controlCases,1)*4,8);
traces = struct([]);
row = 0;
for j = 1:numel(testData)
    ref = testData{j}.ref;
    dist = testData{j}.dist;
    [yPart1,uPart1] = simulate_pdc(P,linearControllers,[],1,ref,dist,ctrl,false);
    [yLinear,uLinear] = simulate_pdc(P,linearControllers,[],1,ref,dist,ctrl,true);
    [yFuzzy,uFuzzy] = simulate_pdc(P,fuzzyControllers,modelTuned.centers, ...
                                   gammaInitial,ref,dist,ctrl,true);
    [yTuned,uTuned] = simulate_pdc(P,fuzzyControllers,controllerCenters, ...
                                   gammaTuned,ref,dist,ctrl,true);

    traces(j).ref = ref;
    traces(j).dist = dist;
    traces(j).y = [yPart1,yLinear,yFuzzy,yTuned];
    traces(j).u = [uPart1,uLinear,uFuzzy,uTuned];
    for method = 1:4
        row = row+1;
        yCase = traces(j).y(:,method);
        uCase = traces(j).u(:,method);
        controlResults(row,:) = [j,method, ...
            mean((yCase(refMask)-ref(refMask)).^2), ...
            mean((yCase(distMask)-ref(distMask)).^2), ...
            mean((yCase(refMask|distMask)-ref(refMask|distMask)).^2), ...
            min(uCase),max(uCase),max(abs(diff([P.u0;uCase])))];
    end
end
plots.pdc(tControl,traces,controlCases);

controlSummary = array2table(controlResults,'VariableNames', ...
    {'Scenariusz','Regulator','MSE_zadanie','MSE_zaklocenie','MSE_laczny', ...
     'Min_u','Max_u','Maks_zmiana_u'});
disp(controlSummary);

totalsPart2 = accumarray(controlSummary.Regulator, ...
    controlSummary.MSE_laczny,[],@mean);
fprintf('Sredni MSE [DMC cz.1, Liniowy ogr., PDC, PDC strojony] = %s cm^2\n', ...
    mat2str(totalsPart2.',5));
fprintf('Regulator PDC: gamma przed = %.2f, po strojeniu = %.2f\n', ...
    gammaInitial,gammaTuned);
fprintf('Wezly regulatora po strojeniu = %s cm; sredni MSE = %.4f cm^2\n', ...
    mat2str(controllerCenters,4),bestControlScore);

controllerWeights = ts_memberships(levels,controllerCenters).^gammaTuned;
controllerWeights = controllerWeights./sum(controllerWeights,2);
plots.plot_control_memberships(levels,initialWeights,controllerWeights,gammaTuned);

writetable(modelSummary,fullfile(outdir,'jakosc_modeli.csv'));
writetable(amplitudeSummary,fullfile(outdir,'jakosc_modeli_amplitudy.csv'));
writetable(controlSummary,fullfile(outdir,'jakosc_regulacji.csv'));
save(fullfile(outdir,'dane_czesc2.mat'),'P','models','modelInitial', ...
     'modelTuned','modelSummary','amplitudeSummary','stepMse','nominalDrift', ...
     'openSteps','truth','validationSteps','validationRmse', ...
     'controlCases','controlSummary','traces','tControl','refStart','distStart', ...
     'gammaCandidates','gammaScores','gammaTuned','controllerCenters','ctrl');

%% Funkcje pomocnicze
function y = simulate_plant_step(P,stepSize,n)
    h = P.h0;
    y = zeros(n+1,1);
    y(1) = h(2);
    for k = 1:n
        du = stepSize(1)*(k>P.delay);
        dd = stepSize(2);
        h = rk4_step(@(state) plant_rhs(P,state,P.u0+du,P.d0+dd),h,P.Ts,5);
        h = max(h,0.05);
        y(k+1) = h(2);
    end
end

function y = simulate_arx_step(P,model,stepSize,n)
    % Symulacja dyskretnego modelu TS (wzor 4.41) w petli otwartej.
    y = zeros(n+1,1);
    y(1) = P.h0(2);
    for k = 1:n
        yk = y(k);
        if k==1, ykm1 = P.h0(2); else, ykm1 = y(k-1); end
        q1 = P.q0 + stepSize(1)*(k>P.delay) + stepSize(2);
        if k==1, q2 = P.q0; else, q2 = P.q0 + stepSize(1)*(k-1>P.delay) + stepSize(2); end
        w = ts_memberships(yk,model.centers);
        b = model.b*w.';
        c = model.c*w.';
        a0q = model.a0q*w.';
        y(k+1) = max(b(1)*yk + b(2)*ykm1 + c(1)*q1 + c(2)*q2 + a0q, 0.05);
    end
end

function dh = plant_rhs(P,h,u,d)
    f2 = P.a1*sqrt(max(h(1),0));
    f3 = P.a2*sqrt(max(h(2),0));
    cap2 = 3*P.C2*max(h(2),0.05)^2;
    dh = [(u+d-f2)/P.A1; (f2-f3)/cap2];
end

function controllers = design_dmc(P,model,ctrl)
    % Lokalne analityczne regulatory DMC (wzory 4.72-4.73) dla kazdej reguly TS.
    S = local_step_responses(P,model,ctrl.D);
    for rule = 1:model.count
        [M,freeMatrix] = dmc_matrices(S(:,rule),ctrl.N,ctrl.Nu);
        K = (M'*M+ctrl.lambda*eye(ctrl.Nu))\M';
        controllers(rule).gain = K(1,:); %#ok<AGROW>
        controllers(rule).freeMatrix = freeMatrix; %#ok<AGROW>
    end
end

function score = pdc_mse(P,controllers,centers,gamma,testData,ctrl,refMask,distMask)
    score = 0;
    allMask = refMask | distMask;
    for j = 1:numel(testData)
        trial = testData{j};
        y = simulate_pdc(P,controllers,centers,gamma,trial.ref,trial.dist,ctrl,true);
        err = y-trial.ref;
        score = score + mean(err(allMask).^2);
    end
    score = score / numel(testData);
end

function [y,u] = simulate_pdc(P,controllers,centers,gamma,ref,dist,ctrl,limited)
    % Analityczny rozmyty DMC (wzory 4.74-4.75) z rzutowaniem na ograniczenia.
    n = numel(ref)-1;
    y = zeros(n+1,1); y(1) = P.h0(2);
    u = zeros(n,1);
    h = P.h0;
    pastDu = zeros(ctrl.D-1,1);
    nRules = numel(controllers);
    for k = 1:n
        if k==1, previousU = P.u0; else, previousU = u(k-1); end
        if nRules==1
            weights = 1;
        else
            weights = ts_memberships(y(k),centers).^gamma;
            weights = weights/sum(weights);
        end
        du = 0;
        for rule = 1:nRules
            free = y(k) + controllers(rule).freeMatrix*pastDu;
            du = du + weights(rule)*(controllers(rule).gain*(ref(k)-free));
        end
        if limited
            lower = max(-ctrl.duMax,ctrl.uMin-previousU);
            upper = min( ctrl.duMax,ctrl.uMax-previousU);
            du = min(max(du,lower),upper);
        end
        u(k) = previousU + du;
        pastDu = [du;pastDu(1:end-1)];

        if k>P.delay, actualU = u(k-P.delay); else, actualU = P.u0; end
        h = rk4_step(@(state) plant_rhs(P,state,actualU,dist(k)),h,P.Ts,5);
        h = max(h,0.05);
        y(k+1) = h(2);
    end
end
