%% Kontrole numeryczne implementacji czesci 3
clear; clc;
load(fullfile('wyniki_czesc3','dane_czesc3.mat'), ...
    'P','modelTuned','ctrl','traces','slTraces','tControl', ...
    'comparisonSummary','horizons','predictionHorizons','lambdaSL','referenceTime','localSteps', ...
    'validationCases','validationTraces','validationSummary','refStart','distStart');

%% Zapisane lokalne odpowiedzi jednostkowe a niezalezna dyskretyzacja ZOH
maxStepError = 0;
Ad = zeros(2,2,modelTuned.count);
Bd = zeros(2,modelTuned.count);
for rule = 1:modelTuned.count
    A = modelTuned.A(:,:,rule);
    E = expm([A,modelTuned.B;zeros(1,3)]*P.Ts);
    Ad(:,:,rule) = E(1:2,1:2);
    Bd(:,rule) = E(1:2,3);
    state = zeros(2,1);
    calculated = zeros(ctrl.D,1);
    for i = 1:numel(calculated)
        if i>P.delay, state=Ad(:,:,rule)*state+Bd(:,rule); end
        calculated(i) = state(2);
    end
    maxStepError = max(maxStepError,max(abs(calculated-localSteps(:,rule))));
end
fprintf('Maks. blad odpowiedzi skokowej wzgledem step: %.3g cm\n',maxStepError);
assert(maxStepError<1e-10);

%% Predykcja przy dwoch aktywnych regulach a symulacja obu modeli stanu
midpoint = mean(modelTuned.centers(2:3));
w = ts_memberships(midpoint,modelTuned.centers);
assert(max(abs(w-[0 .5 .5 0]))<1e-12);
s = localSteps*w.';
Nu = horizons(end);
N = predictionHorizons(end);
[M,F] = dmc_matrices(s,N,Nu);
assert(max(abs(F(:,end)-(s(end)-s(end-1))))<1e-12);
pastDu = zeros(ctrl.D-1,1);
% Wszystkie skoki sa dostatecznie mlode, aby ominac obciecie ogona po D.
pastDu([1,4,17,69,140,185]) = [1;-2;3;1;-1;2];
moves = 2*(-0.6).^(0:Nu-1)';
k = ctrl.D+5;
input = zeros(k+N-1,1);
for j = 1:numel(pastDu), input(k-j:end)=input(k-j:end)+pastDu(j); end
for j = 1:Nu, input(k+j-1:end)=input(k+j-1:end)+moves(j); end
state = zeros(2,modelTuned.count);
output = zeros(numel(input)+1,1);
for m = 1:numel(input)
    if m>P.delay, delayedInput=input(m-P.delay);
    else, delayedInput=0; end
    for rule = 1:modelTuned.count
        state(:,rule) = Ad(:,:,rule)*state(:,rule)+Bd(:,rule)*delayedInput;
    end
    output(m+1)=state(2,:)*w.';
end
predicted = output(k)+F*pastDu+M*moves;
matrixError = max(abs(predicted-output(k+(1:N))));
fprintf('Maks. blad predykcji DMC wzgledem symulacji stanu: %.3g cm\n',matrixError);
assert(matrixError<1e-10);

%% Wagi TS w zakresie modelu i na obu ramionach
grid = linspace(modelTuned.centers(1)-1,modelTuned.centers(end)+1,301)';
weights = ts_memberships(grid,modelTuned.centers);
weightError = max(abs(sum(weights,2)-1));
fprintf('Maks. blad sumy wag TS: %.3g\n',weightError);
assert(weightError<1e-12 && all(weights(:)>=0));

%% Niezalezna ocena wszystkich prob projektowych i walidacyjnych
saved = [slTraces;validationTraces];
refs = {traces.ref}; disturbances = {traces.dist};
for j = 1:size(validationCases,1)
    refs{end+1} = P.h0(2)+validationCases(j,1)*(tControl>=refStart);
    disturbances{end+1} = P.d0+validationCases(j,2)*(tControl>=distStart);
end
task = tControl>=refStart & tControl<distStart;
disturbanceMask = tControl>=distStart;
maxConstraintError = 0; maxMetricError = 0;
moveError = zeros(numel(horizons),1);
lsOptions = optimoptions('lsqlin','Algorithm','active-set','Display','off');
for scenario = 1:size(saved,1)
    ref = refs{scenario};
    for v = 1:numel(horizons)
        y = saved(scenario,v).y; u = saved(scenario,v).u;
        du = diff([P.u0;u]);
        maxConstraintError = max([maxConstraintError;u-ctrl.uMax;ctrl.uMin-u;abs(du)-ctrl.duMax]);
        if scenario<=numel(traces)
            row = comparisonSummary(comparisonSummary.Scenariusz==scenario & comparisonSummary.Regulator==v+2,:);
        else
            row = validationSummary(validationSummary.Scenariusz==scenario-numel(traces) & validationSummary.Wariant==v,:);
        end
        err = y-ref;
        signStep = sign(ref(end)-ref(1));
        directed = signStep*y(task);
        % Pary i<=j daja wszystkie cofniecia od wczesniejszej probki i do pozniejszej j.
        retreat = max(triu(directed-directed.'),[],'all');
        expected = [P.Ts*sum(abs(err(task))),P.Ts*sum(abs(err(disturbanceMask))), ...
            min(u),max(u),max(abs(du)),sum(abs(du)),retreat, ...
            max(0,max(signStep*err(task))),max(abs(err(disturbanceMask)))];
        measured = row{1,{'IAE_zadanie','IAE_zaklocenie','Min_u','Max_u','Maks_zmiana_u', ...
            'Suma_zmian_u','Cofniecie_cm','Przeregulowanie_cm','Maks_blad_po_FD_cm'}};
        maxMetricError = max(maxMetricError,max(abs(expected-measured)));

        % Niezalezna rekonstrukcja celu i rozwiazania QP z zapisanych y i u.
        target = ref;
        if referenceTime(v)>0
            pole = exp(-P.Ts/referenceTime(v));
            target = filter(1-pole,[1,-pole],ref,pole*ref(1));
        end
        Nu = horizons(v);
        L = tril(ones(Nu)); onesU = ones(Nu,1); limits = ctrl.duMax*onesU;
        history = zeros(ctrl.D-1,1); prev = P.u0;
        for k = 1:numel(u)
            w = ts_memberships(y(k),modelTuned.centers);
            [M,F] = dmc_matrices(localSteps*w.',predictionHorizons(v),Nu);
            rhs = target(k)-(y(k)+F*history);
            if Nu==1
                raw = M'*rhs/(M'*M+lambdaSL(v));
                move = min([max([raw,-ctrl.duMax,ctrl.uMin-prev]),ctrl.duMax,ctrl.uMax-prev]);
            else
                G = [M;sqrt(lambdaSL(v))*eye(Nu)]; b = [rhs;zeros(Nu,1)];
                plan = G\b;
                values = prev+L*plan;
                if any(abs(plan)>limits) || any(values<ctrl.uMin | values>ctrl.uMax)
                    bounds = [(ctrl.uMax-prev)*onesU;(prev-ctrl.uMin)*onesU];
                    [plan,~,~,flag] = lsqlin(G,b,[L;-L],bounds,[],[],-limits,limits,zeros(Nu,1),lsOptions);
                    assert(flag>0);
                end
                move = plan(1);
            end
            moveError(v) = max(moveError(v),abs(move-du(k)));
            history = [du(k);history(1:end-1)];
            prev = u(k);
        end
    end
end
fprintf('Najwieksze naruszenie ograniczen: %.3g; blad wskaznikow: %.3g\n',maxConstraintError,maxMetricError);
fprintf('Maks. blad ruchu QP w kolejnych wariantach SL: %s\n',mat2str(moveError.',3));
assert(maxConstraintError<1e-8 && maxMetricError<1e-8 && all(moveError<1e-6));

%% Niezalezny bilans objetosci dla zapisanego sterowania
maxPlantError = 0;
for scenario = 1:size(saved,1)
    for variant = 1:numel(horizons)
        control = saved(scenario,variant).u;
        disturbance = disturbances{scenario};
        volume = [P.A1*P.h0(1);P.C2*P.h0(2)^3];
        checkY = zeros(numel(control)+1,1);
        checkY(1) = P.h0(2);
        for k = 1:numel(control)
            if k>P.delay, delayedU=control(k-P.delay);
            else, delayedU=P.u0; end
            dt = P.Ts/5;
            for substep = 1:5
                a = volume_rhs(P,volume,delayedU,disturbance(k));
                b = volume_rhs(P,volume+dt*a/2,delayedU,disturbance(k));
                c = volume_rhs(P,volume+dt*b/2,delayedU,disturbance(k));
                d = volume_rhs(P,volume+dt*c,delayedU,disturbance(k));
                volume = volume+dt*(a+2*b+2*c+d)/6;
            end
            checkY(k+1) = (volume(2)/P.C2)^(1/3);
        end
        maxPlantError = max(maxPlantError, ...
            max(abs(checkY-saved(scenario,variant).y)));
    end
end
fprintf('Maks. roznica symulacji wysokosci i objetosci: %.3g cm\n',maxPlantError);
assert(maxPlantError<1e-4);

function dv = volume_rhs(P,v,u,d)
    h1 = max(v(1)/P.A1,0);
    h2 = max(v(2)/P.C2,0)^(1/3);
    flow1 = P.a1*sqrt(h1);
    flow2 = P.a2*sqrt(h2);
    dv = [u+d-flow1;flow1-flow2];
end
