%% SZAU, zadanie 1, czesc 2: modele Takagiego-Sugeno i rozmyty DMC
% Uruchom z folderu projektu: run_part2
% Wymaga Control System Toolbox i Fuzzy Logic Toolbox (trimf).
clear; close all; clc;

%% Parametry obiektu (takie same jak w czesci 1)
P = tank_parameters();

outdir = fullfile(fileparts(mfilename('fullpath')),'wyniki_czesc2');
plots = ResultPlots(outdir);

%% Modele TS: od 2 do 5 modeli lokalnych
% h2 jako zmienna reguly: jest mierzona i silnie zmienia pojemnosc V2'(h2).
% Zakres odpowiada F1in = 50..150 przy FD = 20 cm^3/s.
hRange = [(70/P.a2)^2, (170/P.a2)^2];
modelCounts = 2:5;
models = cell(size(modelCounts));
for j = 1:numel(modelCounts)
    centers = linspace(hRange(1),hRange(2),modelCounts(j));
    models{j} = make_ts_model(P,centers);
end
linearModel = make_ts_model(P,P.h0(2));

%% Ocena modeli: skoki F1in i FD z punktu pracy
% Kolumny: Delta F1in, Delta FD. Kazda proba zaczyna sie od rownowagi.
openSteps = [-50 0; -25 0; -10 0; 10 0; 25 0; 50 0; ...
              0 -10; 0 -5; 0 -2; 0 2; 0 5; 0 10];
nOpen = 300;
tOpen = (0:nOpen)'*P.Ts;
truth = zeros(nOpen+1,size(openSteps,1));
for j = 1:size(openSteps,1)
    truth(:,j) = simulate_step(P,[],openSteps(j,:),nOpen);
end

modelRmse = zeros(numel(models),1);
modelMaxError = zeros(numel(models),1);
% Kolumny: liniowy, TS2, TS3, TS4, TS5 i wybrany TS4 po strojeniu.
stepMse = zeros(size(openSteps,1),numel(models)+2);
modelOutputs = zeros([size(truth),numel(models)]);
linearOutputs = zeros(size(truth));
linearErrors = zeros(size(truth));
for j = 1:size(openSteps,1)
    linearOutputs(:,j) = simulate_step(P,linearModel,openSteps(j,:),nOpen);
    linearErrors(:,j) = linearOutputs(:,j)-truth(:,j);
end
linearRmse = sqrt(mean(linearErrors(:).^2));
linearMaxError = max(abs(linearErrors(:)));
stepMse(:,1) = mean(linearErrors.^2,1).';
for m = 1:numel(models)
    errors = zeros(size(truth));
    for j = 1:size(openSteps,1)
        yTS = simulate_step(P,models{m},openSteps(j,:),nOpen);
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

% Wybrany model do dalszych badan: 4 reguly (dobry kompromis jakosci i liczby regul).
selectedCount = 4;
selectedIndex = find(modelCounts==selectedCount);
modelInitial = models{selectedIndex};

% Strojenie: przesuwamy dwa wewnetrzne wezly, pozostawiajac konce zakresu.
middle = modelInitial.centers(2:3);
shifts = [-3 -1.5 0 1.5 3];
bestRmse = inf;
modelTuned = modelInitial;
tunedOutputs = modelOutputs(:,:,selectedIndex);
candidateOutputs = zeros(size(truth));
for a = shifts
    for b = shifts
        candidateCenters = [hRange(1),middle(1)+a,middle(2)+b,hRange(2)];
        if any(diff(candidateCenters)<=0), continue; end
        % Macierze lokalne pozostaja stale; zmieniaja sie tylko wagi regul.
        candidateModel = modelInitial;
        candidateModel.centers = candidateCenters;
        candidateMse = zeros(size(openSteps,1),1);
        for j = 1:size(openSteps,1)
            yTS = simulate_step(P,candidateModel,openSteps(j,:),nOpen);
            error = yTS-truth(:,j);
            candidateOutputs(:,j) = yTS;
            candidateMse(j) = mean(error.^2);
        end
        rmse = sqrt(mean(candidateMse));
        % Wybierz najlepszy model dla zestawu 12 skokow projektowych.
        if rmse<bestRmse
            bestRmse = rmse;
            modelTuned = candidateModel;
            tunedOutputs = candidateOutputs;
            stepMse(:,end) = candidateMse;
        end
    end
end
if isinf(bestRmse), bestRmse=modelRmse(selectedIndex); end
fprintf('Model TS 4: RMSE przed = %.4f cm, po = %.4f cm; wezly = %s cm\n', ...
    modelRmse(selectedIndex),bestRmse,mat2str(modelTuned.centers,4));

% Osobny zestaw skokow nieuzyty przy wyborze wezlow.
validationSteps = [-40 0;40 0;0 -7;0 7];
validationError = zeros(nOpen+1,size(validationSteps,1),2);
for j = 1:size(validationSteps,1)
    yTrue = simulate_step(P,[],validationSteps(j,:),nOpen);
    validationError(:,j,1) = simulate_step(P,modelInitial, ...
        validationSteps(j,:),nOpen)-yTrue;
    validationError(:,j,2) = simulate_step(P,modelTuned, ...
        validationSteps(j,:),nOpen)-yTrue;
end
validationRmse = [sqrt(mean(validationError(:,:,1).^2,'all')); ...
                  sqrt(mean(validationError(:,:,2).^2,'all'))];
fprintf('Dodatkowe skoki TS 4: RMSE przed = %.4f cm, po = %.4f cm\n', ...
    validationRmse(1),validationRmse(2));

% W kazdej grupie sa cztery rownej dlugosci odpowiedzi: oba wejscia i oba znaki.
amplitudePct = 100*max(abs(openSteps)./[P.u0,P.d0],[],2);
levels = unique(amplitudePct);
groupRmse = zeros(numel(levels),size(stepMse,2));
for j = 1:numel(levels)
    groupRmse(j,:) = sqrt(mean(stepMse(amplitudePct==levels(j),:),1));
end
amplitudeSummary = array2table([levels,groupRmse],'VariableNames', ...
    {'Skok_proc','Liniowy_cm','TS2_cm','TS3_cm','TS4_cm','TS5_cm','TS4_strojony_cm'});
disp(amplitudeSummary);

% Zmiana poziomu po jednym opoznieniu transportowym przy stalych obu doplywach.
driftModels = [{linearModel},models,{modelTuned}];
nominalDrift = zeros(1,numel(driftModels));
for j = 1:numel(driftModels)
    yNominal = simulate_step(P,driftModels{j},[0,0],P.delay);
    nominalDrift(j) = yNominal(end)-P.h0(2);
end

%% Wykresy modeli i funkcji przynaleznosci
levels = linspace(hRange(1),hRange(2),500)';
modelWeights = cellfun(@(m) ts_memberships(levels,m.centers),models,'UniformOutput',false);
initialWeights = ts_memberships(levels,modelTuned.centers);
plots.plot_memberships(levels,modelWeights,modelCounts,initialWeights);
plots.plot_model_responses(modelOutputs(:,:,selectedIndex),tunedOutputs,openSteps,truth,tOpen);
plots.plot_all_models(linearOutputs,modelOutputs,modelCounts,openSteps,truth,tOpen);

%% Analityczny rozmyty DMC z ograniczeniami
% Nu=1 pozwala wyznaczyc kazdy lokalny ruch analitycznie i ograniczyc go
% dokladnie do przedzialu dopuszczalnego. Wynik PDC to suma wazona ruchow.
ctrl.N = 50;
ctrl.D = 250;
ctrl.lambda = 0.05;
ctrl.uMin = 0;         % cm^3/s
ctrl.uMax = 200;       % cm^3/s
ctrl.duMax = 30;       % cm^3/s na probke

linearConstrained = design_dmc(P,linearModel,ctrl,1);
fuzzyDmc = design_dmc(P,modelTuned,ctrl,1);
linearPart1 = design_dmc(P,linearModel,ctrl,10);

% Scenariusze: zmiana zadania w 400 s, zmiana FD w 3200 s.
controlCases = [ 2,  10; ...
                 4,  20; ...
                 8,  30; ...
                -8, -10];
tControl = (0:320)'*P.Ts;
testData = cell(size(controlCases,1),1);
for j = 1:numel(testData)
    testData{j}.ref = P.h0(2)+controlCases(j,1)*(tControl>=400);
    testData{j}.dist = P.d0+controlCases(j,2)*(tControl>=3200);
end

% PDC: regulator uzywa tych samych funkcji przynaleznosci co wybrany model.
gammaInitial = 1;
gammaCandidates = [0.5 0.75 1 1.25 1.5 2 3];
gammaScores = zeros(size(gammaCandidates));
for g = 1:numel(gammaCandidates)
    gammaScores(g) = pdc_iae(P,fuzzyDmc,modelTuned.centers,gammaCandidates(g),testData,ctrl);
end
[~,bestGammaIndex] = min(gammaScores);
gammaTuned = gammaCandidates(bestGammaIndex);
controllerCenters = modelTuned.centers;
bestControlScore = gammaScores(bestGammaIndex);
centerShifts = [-2 0 2];
for firstShift = centerShifts
    for secondShift = centerShifts
        candidateCenters = modelTuned.centers;
        candidateCenters(2:3) = candidateCenters(2:3) ...
                              + [firstShift,secondShift];
        if any(diff(candidateCenters)<=0), continue; end
        score = pdc_iae(P,fuzzyDmc,candidateCenters,gammaTuned,testData,ctrl);
        if score<bestControlScore
            bestControlScore = score;
            controllerCenters = candidateCenters;
        end
    end
end

% Kolumny: scenariusz, regulator, IAE zadania, IAE zaklocenia, min/max u.
controlResults = zeros(size(controlCases,1)*4,7);
row = 0;

for j = 1:numel(testData)
    ref = testData{j}.ref;
    dist = testData{j}.dist;
    [yPart1,uPart1] = simulate_control(P,linearPart1,[],1,ref,dist,ctrl,false);
    [yLinear,uLinear] = simulate_control(P,linearConstrained,[],1,ref,dist,ctrl,true);
    [yFuzzy,uFuzzy] = simulate_control(P,fuzzyDmc,modelTuned.centers, ...
                                     gammaInitial,ref,dist,ctrl,true);
    [yTuned,uTuned] = simulate_control(P,fuzzyDmc,controllerCenters, ...
                                     gammaTuned,ref,dist,ctrl,true);

    traces(j).ref = ref; traces(j).dist = dist;
    traces(j).y = [yPart1,yLinear,yFuzzy,yTuned];
    traces(j).u = [uPart1,uLinear,uFuzzy,uTuned];
    for method = 1:4
        row = row+1;
        yCase = traces(j).y(:,method);
        uCase = traces(j).u(:,method);
        refMask = tControl>=400 & tControl<3200;
        distMask = tControl>=3200;
        controlResults(row,:) = [j,method, ...
            P.Ts*sum(abs(yCase(refMask)-ref(refMask))), ...
            P.Ts*sum(abs(yCase(distMask)-ref(distMask))), ...
            min(uCase),max(uCase),max(abs(diff([P.u0;uCase])))];
    end

end
plots.pdc(tControl,traces,controlCases);

controlSummary = array2table(controlResults,'VariableNames', ...
    {'Scenariusz','Regulator','IAE_zadanie','IAE_zaklocenie', ...
     'Min_u','Max_u','Maks_zmiana_u'});
disp(controlSummary);
fprintf('Regulator: gamma przed = %.2f, po strojeniu = %.2f\n', ...
        gammaInitial,gammaTuned);
fprintf('Suma IAE dla gamma = %s\n',mat2str(gammaScores,5));
fprintf('Wezly regulatora po strojeniu = %s cm; suma IAE = %.1f cm*s\n', ...
        mat2str(controllerCenters,4),bestControlScore);

controllerWeights = ts_memberships(levels,controllerCenters).^gammaTuned;
plots.plot_control_memberships(levels,initialWeights,controllerWeights./sum(controllerWeights,2),gammaTuned);
writetable(modelSummary,fullfile(outdir,'jakosc_modeli.csv'));
writetable(amplitudeSummary,fullfile(outdir,'jakosc_modeli_amplitudy.csv'));
writetable(controlSummary,fullfile(outdir,'jakosc_regulacji.csv'));
save(fullfile(outdir,'dane_czesc2.mat'),'P','models','modelInitial', ...
     'modelTuned','modelSummary','amplitudeSummary','stepMse','nominalDrift', ...
     'openSteps','truth','validationSteps', ...
     'validationRmse','controlCases', ...
     'controlSummary','traces','tControl','gammaCandidates','gammaScores', ...
     'gammaTuned','controllerCenters','ctrl');

%% Funkcje pomocnicze
function model = make_ts_model(P,centers)
    model.centers = centers(:)';
    model.count = numel(centers);
    model.A = zeros(2,2,model.count);
    model.B = [1/P.A1;0];
    model.heq = zeros(2,model.count);
    model.qeq = zeros(1,model.count);
    for i = 1:model.count
        q = P.a2*sqrt(centers(i));
        h1 = (q/P.a1)^2;
        k1 = P.a1/(2*sqrt(h1));
        k2 = P.a2/(2*sqrt(centers(i)));
        cap2 = 3*P.C2*centers(i)^2;
        model.heq(:,i) = [h1;centers(i)];
        model.qeq(i) = q;
        model.A(:,:,i) = [-k1/P.A1,0; k1/cap2,-k2/cap2];
    end
end

function y = simulate_step(P,model,stepSize,n)
    h = P.h0;
    y = zeros(n+1,1);
    y(1) = h(2);
    for k = 1:n
        % Skok FD dziala od razu, a skok F1in po opoznieniu tau.
        du = stepSize(1)*(k>P.delay);
        dd = stepSize(2);
        h = rk4_step(@(state) model_rhs(P,model,state,du,dd),h,P.Ts,5);
        h = max(h,0.05);
        y(k+1) = h(2);
    end
end

function dh = model_rhs(P,model,h,du,dd)
    if isempty(model)
        f2 = P.a1*sqrt(max(h(1),0));
        f3 = P.a2*sqrt(max(h(2),0));
        cap2 = 3*P.C2*max(h(2),0.05)^2;
        dh = [(P.q0+du+dd-f2)/P.A1;(f2-f3)/cap2];
    else
        weights = ts_memberships(h(2),model.centers);
        dh = zeros(2,1);
        for i = 1:model.count
            local = model.A(:,:,i)*(h-model.heq(:,i)) ...
                  +model.B*(P.q0+du+dd-model.qeq(i));
            dh = dh+weights(i)*local;
        end
    end
end

function controllers = design_dmc(P,model,ctrl,Nu)
    % Jeden lokalny analityczny DMC na kazda regule modelu TS.
    S = local_step_responses(P,model,ctrl.D);
    for rule = 1:model.count
        [M,freeMatrix] = dmc_matrices(S(:,rule),ctrl.N,Nu);
        gain = (M'*M+ctrl.lambda*eye(Nu))\M';
        controllers(rule).gain = gain(1,:);
        controllers(rule).freeMatrix = freeMatrix;
    end
end

function score = pdc_iae(P,controllers,centers,gamma,testData,ctrl)
    score = 0;
    for j = 1:numel(testData)
        trial = testData{j};
        y = simulate_control(P,controllers,centers,gamma,trial.ref,trial.dist,ctrl,true);
        score = score+P.Ts*sum(abs(y-trial.ref));
    end
end

function [y,u] = simulate_control(P,controllers,centers,gamma,ref,dist,ctrl,limited)
    n = numel(ref)-1;
    y = zeros(n+1,1); y(1) = P.h0(2);
    u = zeros(n,1);
    h = P.h0;
    pastDu = zeros(ctrl.D-1,1);
    for k = 1:n
        if k==1, previousU=P.u0; else, previousU=u(k-1); end
        if numel(controllers)==1
            weights = 1;
        else
            weights = ts_memberships(y(k),centers).^gamma;
            weights = weights/sum(weights);
        end

        localMove = zeros(1,numel(controllers));
        for rule = 1:numel(controllers)
            free = y(k)+controllers(rule).freeMatrix*pastDu;
            du = controllers(rule).gain*(ref(k)-free);
            if limited
                lower = max(-ctrl.duMax,ctrl.uMin-previousU);
                upper = min(ctrl.duMax,ctrl.uMax-previousU);
                du = min(max(du,lower),upper);
            end
            localMove(rule) = du;
        end
        du = sum(weights.*localMove);
        u(k) = previousU+du;
        pastDu = [du;pastDu(1:end-1)];

        if k>P.delay, actualU=u(k-P.delay); else, actualU=P.u0; end
        inputChange = [actualU-P.u0,dist(k)-P.d0];
        h = rk4_step(@(state) model_rhs(P,[],state, ...
            inputChange(1),inputChange(2)),h,P.Ts,5);
        h = max(h,0.05);
        y(k+1) = h(2);
    end
end
