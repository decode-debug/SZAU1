%% SZAU, zadanie 1, czesc 2: dyskretne modele Takagiego-Sugeno i rozmyty analityczny DMC (PDC)
% Uruchom z folderu projektu: run_part2
% Wymaga Control System Toolbox (ss, c2d, step, tfdata) i Optimization Toolbox (fmincon).
clear; close all; clc;

%% Parametry obiektu (identyczne jak w czesci 1)
P = tank_parameters();
outdir = fullfile(fileparts(mfilename('fullpath')),'wyniki_czesc2');
plots = ResultPlots(outdir);

%% Modele rozmyte TS: od 2 do 5 modeli lokalnych
% Zmienna przeslanek: mierzony poziom h2 (zakres odpowiada doplywowi q = 50..180 cm^3/s).
hRange = [(50/P.a2)^2, (180/P.a2)^2];
modelCounts = 2:5;
models = arrayfun(@(c) build_arx_model(P,linspace(hRange(1),hRange(2),c)), ...
    modelCounts,'UniformOutput',false);
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

stepMse = zeros(size(openSteps,1),numel(models)+2);
[linearOutputs,stepMse(:,1),linearRmse,linearMaxError] = ...
    eval_arx_steps(P,linearModel,openSteps,nOpen,truth);

modelOutputs = zeros([size(truth),numel(models)]);
modelRmse = zeros(numel(models),1);
modelMaxError = zeros(numel(models),1);
for m = 1:numel(models)
    [modelOutputs(:,:,m),stepMse(:,m+1),modelRmse(m),modelMaxError(m)] = ...
        eval_arx_steps(P,models{m},openSteps,nOpen,truth);
end
modelSummary = table([1;modelCounts'],[linearRmse;modelRmse], ...
    [linearMaxError;modelMaxError], ...
    'VariableNames',{'Liczba_regul','RMSE_cm','Maks_blad_cm'});
disp(modelSummary);

%% Wybor modelu (5 regul) i strojenie wezlow funkcji przynaleznosci
selectedIndex = find(modelCounts==5);
modelInitial = models{selectedIndex};
modelTuned = modelInitial;

% Siatka przesuniec trzech wezlow wewnetrznych, a nastepnie dostrojenie fmincon.
middle = modelInitial.centers(2:4);
shifts = [-2 -1 0 1 2];
bestMse = inf;
for s1 = shifts
    for s2 = shifts
        for s3 = shifts
            candCenters = [hRange(1),middle+[s1,s2,s3],hRange(2)];
            if any(diff(candCenters)<=2), continue; end
            mseVal = model_ts_mse(P,modelInitial,candCenters,openSteps,nOpen,truth);
            if mseVal<bestMse
                bestMse = mseVal;
                modelTuned.centers = candCenters;
            end
        end
    end
end

objModel = @(cIn) model_ts_mse(P,modelInitial,[hRange(1),cIn(:)',hRange(2)],openSteps,nOpen,truth);
optSq = optimoptions('fmincon','Display','off','Algorithm','sqp');
bestInner = fmincon(objModel,modelTuned.centers(2:4),[1 -1 0; 0 1 -1],[-1.5; -1.5], ...
    [],[],(hRange(1)+1.5)*ones(1,3),(hRange(2)-1.5)*ones(1,3),[],optSq);
modelTuned.centers = [hRange(1),bestInner(:)',hRange(2)];

[tunedOutputs,stepMse(:,end),bestRmse] = eval_arx_steps(P,modelTuned,openSteps,nOpen,truth);
fprintf('Model TS 5: MSE przed = %.4f cm^2, po = %.4f cm^2; wezly = %s cm\n', ...
    modelRmse(selectedIndex)^2,bestRmse^2,mat2str(modelTuned.centers,4));

% Walidacja modelu na 4 dodatkowych skokach.
validationSteps = [-40 0; 40 0; 0 -7; 0 7];
valTruth = zeros(nOpen+1,size(validationSteps,1));
for j = 1:size(validationSteps,1)
    valTruth(:,j) = simulate_plant_step(P,validationSteps(j,:),nOpen);
end
[~,~,valBefore] = eval_arx_steps(P,modelInitial,validationSteps,nOpen,valTruth);
[~,~,valAfter]  = eval_arx_steps(P,modelTuned,validationSteps,nOpen,valTruth);
validationRmse = [valBefore; valAfter];
fprintf('Walidacja TS 5: RMSE przed = %.4f cm, po = %.4f cm\n',valBefore,valAfter);

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

%% Wykresy modeli i funkcji przynaleznosci
levels = linspace(hRange(1),hRange(2),500)';
modelWeights = cellfun(@(m) ts_memberships(levels,m.centers),models,'UniformOutput',false);
initialWeights = ts_memberships(levels,modelTuned.centers);
plots.plot_memberships(levels,modelWeights,modelCounts,initialWeights);
plots.plot_model_responses(modelOutputs(:,:,selectedIndex),tunedOutputs,openSteps,truth,tOpen);
plots.plot_all_models(linearOutputs,modelOutputs,modelCounts,openSteps,truth,tOpen);

%% Analityczny rozmyty DMC (PDC) z ograniczeniami sterowania
ctrl.D = 200;
ctrl.N = 50;
ctrl.Nu = 10;
ctrl.lambda = 0.05;
ctrl.uMin = 0;         % cm^3/s
ctrl.uMax = 200;       % cm^3/s
ctrl.duMax = 15;       % cm^3/s na probke

linearControllers = design_dmc(P,linearModel,ctrl);
fuzzyControllers = design_dmc(P,modelTuned,ctrl);

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

% Strojenie funkcji przynaleznosci regulatora PDC (dzwonowe Gaussa o wspolnym sigma).
sigmaCandidates = [2.0 2.5 3.0 3.5 4.0 4.5];
sigmaScores = arrayfun(@(s) pdc_mse(P,fuzzyControllers,modelTuned.centers,s,testData,ctrl,refMask,distMask), ...
    sigmaCandidates);
[bestControlScore,bestSigmaIndex] = min(sigmaScores);
sigmaTuned = sigmaCandidates(bestSigmaIndex);
controllerCenters = modelTuned.centers;

% Siatka przesuniec wezlow wewnetrznych, nastepnie optymalizacja ciagla fmincon.
cShifts = [-2 0 2];
for d1 = cShifts
    for d2 = cShifts
        for d3 = cShifts
            candCenters = modelTuned.centers;
            candCenters(2:4) = candCenters(2:4) + [d1,d2,d3];
            if any(diff(candCenters)<4.0), continue; end
            score = pdc_mse(P,fuzzyControllers,candCenters,sigmaTuned,testData,ctrl,refMask,distMask);
            if score<bestControlScore
                bestControlScore = score;
                controllerCenters = candCenters;
            end
        end
    end
end

p0 = [controllerCenters(2:4), sigmaTuned];
lb = [controllerCenters(1)+4.0, controllerCenters(1)+8.0, controllerCenters(1)+12.0, 2.2];
ub = [controllerCenters(end)-12.0, controllerCenters(end)-8.0, controllerCenters(end)-4.0, 4.5];
objCtrl = @(p) pdc_mse(P,fuzzyControllers,[controllerCenters(1), p(1:3), controllerCenters(end)], ...
    p(4),testData,ctrl,refMask,distMask);
optsCtrl = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxIterations',60,'StepTolerance',1e-4,'OptimalityTolerance',1e-5);
[pOpt,scoreOpt] = fmincon(objCtrl,p0,[1 -1 0 0; 0 1 -1 0],[-4.0; -4.0],[],[],lb,ub,[],optsCtrl);
if scoreOpt<bestControlScore
    bestControlScore = scoreOpt;
    controllerCenters(2:4) = pOpt(1:3);
    sigmaTuned = pOpt(4);
end

%% Symulacje porownawcze 4 regulatorow czesci 2
controlResults = zeros(size(controlCases,1)*4,8);
traces = struct([]);
row = 0;
for j = 1:numel(testData)
    ref = testData{j}.ref;
    dist = testData{j}.dist;
    [yPart1,uPart1]   = simulate_pdc(P,linearControllers,[],0,ref,dist,ctrl,false);
    [yLinear,uLinear] = simulate_pdc(P,linearControllers,[],0,ref,dist,ctrl,true);
    [yFuzzy,uFuzzy]   = simulate_pdc(P,fuzzyControllers,modelTuned.centers,0,ref,dist,ctrl,true);
    [yTuned,uTuned]   = simulate_pdc(P,fuzzyControllers,controllerCenters,sigmaTuned,ref,dist,ctrl,true);

    traces(j).ref = ref;
    traces(j).dist = dist;
    traces(j).y = [yPart1,yLinear,yFuzzy,yTuned];
    traces(j).u = [uPart1,uLinear,uFuzzy,uTuned];
    for method = 1:4
        row = row+1;
        m = control_metrics(P,traces(j).y(:,method),traces(j).u(:,method),ref,refMask,distMask);
        controlResults(row,:) = [j,method,m(1:3),m(6:8)];
    end
end
plots.pdc(tControl,traces,controlCases,refStart,distStart);

controlSummary = array2table(controlResults,'VariableNames', ...
    {'Scenariusz','Regulator','MSE_zadanie','MSE_zaklocenie','MSE_laczny', ...
     'Min_u','Max_u','Maks_zmiana_u'});
disp(controlSummary);

totalsPart2 = accumarray(controlSummary.Regulator,controlSummary.MSE_laczny,[],@mean);
fprintf('Sredni MSE [DMC cz.1, Liniowy ogr., PDC, PDC strojony] = %s cm^2\n', ...
    mat2str(totalsPart2.',5));
fprintf('Regulator PDC po strojeniu (dzwonowe Gaussa): sigma = %.2f cm\n',sigmaTuned);
fprintf('Wezly regulatora po strojeniu = %s cm; sredni MSE = %.4f cm^2\n', ...
    mat2str(controllerCenters,4),bestControlScore);

[~,controllerMu] = ts_memberships(levels,controllerCenters,sigmaTuned);
plots.plot_control_memberships(levels,initialWeights,controllerMu,sigmaTuned);

writetable(modelSummary,fullfile(outdir,'jakosc_modeli.csv'));
writetable(amplitudeSummary,fullfile(outdir,'jakosc_modeli_amplitudy.csv'));
writetable(controlSummary,fullfile(outdir,'jakosc_regulacji.csv'));
save(fullfile(outdir,'dane_czesc2.mat'),'P','models','modelInitial', ...
     'modelTuned','modelSummary','amplitudeSummary','stepMse', ...
     'openSteps','truth','validationSteps','validationRmse', ...
     'controlCases','controlSummary','traces','tControl','refStart','distStart', ...
     'sigmaCandidates','sigmaScores','sigmaTuned','controllerCenters','ctrl');

%% Funkcje pomocnicze
function model = build_arx_model(P,centers)
    % Lokalne linearyzacje i dyskretyzacja ZOH w postaci rownan roznicowych TS.
    model.centers = centers(:)';
    model.count = numel(centers);
    model.A = zeros(2,2,model.count);
    model.B = [1/P.A1;0];
    model.b = zeros(2,model.count);
    model.c = zeros(2,model.count);
    model.a0q = zeros(1,model.count);
    model.a0 = zeros(1,model.count);

    for i = 1:model.count
        h2 = centers(i);
        q = P.a2*sqrt(h2);
        h1 = (q/P.a1)^2;
        k1 = P.a1/(2*sqrt(h1));
        k2 = P.a2/(2*sqrt(h2));
        cap2 = 3*P.C2*h2^2;
        A = [-k1/P.A1, 0; k1/cap2, -k2/cap2];
        model.A(:,:,i) = A;

        sysd = c2d(ss(A,model.B,[0 1],0),P.Ts,'zoh');
        [num,den] = tfdata(sysd,'v');
        b1 = -den(2); b2 = -den(3);
        c1 = num(2);  c2 = num(3);
        model.b(:,i) = [b1;b2];
        model.c(:,i) = [c1;c2];
        model.a0q(i) = (1-b1-b2)*h2 - (c1+c2)*q;
        model.a0(i)  = (1-b1-b2)*h2 - (c1+c2)*(q-P.d0);
    end
end

function [outputs,msePerStep,rmseAll,maxErr] = eval_arx_steps(P,model,openSteps,nOpen,truth)
    nSteps = size(openSteps,1);
    outputs = zeros(nOpen+1,nSteps);
    errors = zeros(nOpen+1,nSteps);
    for j = 1:nSteps
        outputs(:,j) = simulate_arx_step(P,model,openSteps(j,:),nOpen);
        errors(:,j) = outputs(:,j) - truth(:,j);
    end
    msePerStep = mean(errors.^2,1).';
    rmseAll = sqrt(mean(errors(:).^2));
    maxErr = max(abs(errors(:)));
end

function mseVal = model_ts_mse(P,baseModel,centers,openSteps,nOpen,truth)
    cand = baseModel;
    cand.centers = centers;
    [~,stepMse] = eval_arx_steps(P,cand,openSteps,nOpen,truth);
    mseVal = mean(stepMse);
end

function y = simulate_plant_step(P,stepSize,n)
    h = P.h0;
    y = zeros(n+1,1);
    y(1) = h(2);
    for k = 1:n
        du = stepSize(1)*(k>P.delay);
        dd = stepSize(2);
        h = plant_step(P,h,P.u0+du,P.d0+dd);
        y(k+1) = h(2);
    end
end

function y = simulate_arx_step(P,model,stepSize,n)
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

function controllers = design_dmc(P,model,ctrl)
    % Lokalne analityczne regulatory DMC dla kazdej reguly TS.
    tStep = (0:ctrl.D)'*P.Ts;
    for rule = 1:model.count
        sys = ss(model.A(:,:,rule),model.B,[0 1],0,'InputDelay',P.tau);
        resp = step(c2d(sys,P.Ts,'zoh'),tStep);
        [M,freeMatrix] = dmc_matrices(resp(2:end),ctrl.N,ctrl.Nu);
        K = (M'*M+ctrl.lambda*eye(ctrl.Nu))\M';
        controllers(rule).gain = K(1,:); %#ok<AGROW>
        controllers(rule).freeMatrix = freeMatrix; %#ok<AGROW>
    end
end

function score = pdc_mse(P,controllers,centers,sigma,testData,ctrl,refMask,distMask)
    score = 0;
    allMask = refMask | distMask;
    for j = 1:numel(testData)
        trial = testData{j};
        y = simulate_pdc(P,controllers,centers,sigma,trial.ref,trial.dist,ctrl,true);
        score = score + mean((y(allMask)-trial.ref(allMask)).^2);
    end
    score = score / numel(testData);
end

function [y,u] = simulate_pdc(P,controllers,centers,sigma,ref,dist,ctrl,limited)
    % Analityczny rozmyty DMC (PDC) z rzutowaniem na ograniczenia.
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
            weights = ts_memberships(y(k),centers,sigma);
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
        h = plant_step(P,h,actualU,dist(k));
        y(k+1) = h(2);
    end
end
