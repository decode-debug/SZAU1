%% Niezalezne kontrole numeryczne implementacji czesci 3 (FGPC-SL)
clear; clc;
load(fullfile('wyniki_czesc3','dane_czesc3.mat'), ...
    'P','modelTuned','ctrl','traces','slTraces','tControl', ...
    'comparisonSummary','N_sl','Nu_sl','lambdaSL', ...
    'validationCases','validationTraces','validationSummary','refStart','distStart');

%% 1. Weryfikacja dyskretnych modeli lokalnych ARX wzgledem niezaleznej dyskretyzacji ZOH
maxStepError = 0;
for rule = 1:modelTuned.count
    A = modelTuned.A(:,:,rule);
    E = expm([A,modelTuned.B;zeros(1,3)]*P.Ts);
    Ad = E(1:2,1:2);
    Bd = E(1:2,3);
    state = zeros(2,1);
    sSS = zeros(N_sl,1);
    sARX = zeros(N_sl,1);
    b = modelTuned.b(:,rule);
    c = modelTuned.c(:,rule);
    for p = 1:N_sl
        if p>P.delay, state = Ad*state + Bd; end
        sSS(p) = state(2);
        y1 = 0; y2 = 0;
        if p>1, y1 = sARX(p-1); end
        if p>2, y2 = sARX(p-2); end
        sARX(p) = b(1)*y1 + b(2)*y2 + c(1)*(p>P.delay) + c(2)*(p>P.delay+1);
    end
    maxStepError = max(maxStepError,max(abs(sSS-sARX)));
end
fprintf('Maks. blad odpowiedzi skokowej ARX wzgledem ZOH state-space: %.3g cm\n',maxStepError);
assert(maxStepError<1e-10);

%% 2. Weryfikacja wag TS (zgodnosc z trimf oraz partycja jednosci)
grid = linspace(modelTuned.centers(1)-2,modelTuned.centers(end)+2,401)';
weights = ts_memberships(grid,modelTuned.centers);
weightError = max(abs(sum(weights,2)-1));
fprintf('Maks. blad sumy wag TS: %.3g\n',weightError);
assert(weightError<1e-12 && all(weights(:)>=0));

%% 3. Niezalezna rekonstrukcja FGPC-SL (lsqlin) oraz weryfikacja ograniczen i wskaznikow
saved = [slTraces;validationTraces];
refs = {traces.ref};
disturbances = {traces.dist};
for j = 1:size(validationCases,1)
    refs{end+1} = P.h0(2)+validationCases(j,1)*(tControl>=refStart); %#ok<SAGROW>
    disturbances{end+1} = P.d0+validationCases(j,2)*(tControl>=distStart); %#ok<SAGROW>
end
taskMask = tControl>=refStart & tControl<distStart;
distMask = tControl>=distStart;

maxConstraintError = 0;
maxMetricError = 0;
moveError = zeros(numel(lambdaSL),1);
lsOptions = optimoptions('lsqlin','Algorithm','active-set','Display','off');
L = tril(ones(Nu_sl));
onesU = ones(Nu_sl,1);
limits = ctrl.duMax*onesU;

for scenario = 1:size(saved,1)
    ref = refs{scenario};
    for v = 1:numel(lambdaSL)
        y = saved(scenario,v).y;
        u = saved(scenario,v).u;
        du = diff([P.u0;u]);
        maxConstraintError = max([maxConstraintError;u-ctrl.uMax;ctrl.uMin-u;abs(du)-ctrl.duMax]);

        if scenario<=numel(traces)
            row = comparisonSummary(comparisonSummary.Scenariusz==scenario & comparisonSummary.Regulator==v+2,:);
        else
            row = validationSummary(validationSummary.Scenariusz==scenario-numel(traces) & validationSummary.Wariant==v,:);
        end
        err = y-ref;
        signStep = sign(ref(end)-ref(1));
        directed = signStep*y(taskMask);
        retreat = max(triu(directed-directed.'),[],'all');
        expected = [P.Ts*sum(abs(err(taskMask))),P.Ts*sum(abs(err(distMask))), ...
            min(u),max(u),max(abs(du)),sum(abs(du)),retreat, ...
            max(0,max(signStep*err(taskMask))),max(abs(err(distMask)))];
        measured = row{1,{'IAE_zadanie','IAE_zaklocenie','Min_u','Max_u','Maks_zmiana_u', ...
            'Suma_zmian_u','Cofniecie_cm','Przeregulowanie_cm','Maks_blad_po_FD_cm'}};
        maxMetricError = max(maxMetricError,max(abs(expected-measured)));

        % Niezalezna rekonstrukcja kroku FGPC-SL (wzor 4.107 + lsqlin)
        for k = 1:numel(u)
            if k==1, prev = P.u0; else, prev = u(k-1); end
            w = ts_memberships(y(k),modelTuned.centers);
            b = modelTuned.b*w.'; c = modelTuned.c*w.'; a0 = modelTuned.a0*w.';
            s = zeros(N_sl,1);
            for p = 1:N_sl
                y1 = 0; y2 = 0;
                if p>1, y1 = s(p-1); end
                if p>2, y2 = s(p-2); end
                s(p) = b(1)*y1 + b(2)*y2 + c(1)*(p>P.delay) + c(2)*(p>P.delay+1);
            end
            M = toeplitz(s,[s(1),zeros(1,Nu_sl-1)]);
            ykm1 = y(max(1,k-1)); ykm2 = y(max(1,k-2));
            dk = y(k) - (b(1)*ykm1 + b(2)*ykm2 + ...
                c(1)*get_u(u,k-1-P.delay,k-1,prev,P.u0) + ...
                c(2)*get_u(u,k-2-P.delay,k-1,prev,P.u0) + a0);
            free = zeros(N_sl,1); y1 = y(k); y2 = ykm1;
            for p = 1:N_sl
                yN = b(1)*y1 + b(2)*y2 + ...
                    c(1)*get_u(u,k+p-1-P.delay,k-1,prev,P.u0) + ...
                    c(2)*get_u(u,k+p-2-P.delay,k-1,prev,P.u0) + a0 + dk;
                free(p) = yN; y2 = y1; y1 = yN;
            end
            G = [M;sqrt(lambdaSL(v))*eye(Nu_sl)];
            rhs = [ref(k)-free;zeros(Nu_sl,1)];
            plan = G\rhs;
            vals = prev + L*plan;
            if any(abs(plan)>limits) || any(vals<ctrl.uMin | vals>ctrl.uMax)
                bounds = [(ctrl.uMax-prev)*onesU;(prev-ctrl.uMin)*onesU];
                [plan,~,~,flag] = lsqlin(G,rhs,[L;-L],bounds,[],[],-limits,limits,zeros(Nu_sl,1),lsOptions);
                assert(flag>0);
            end
            moveError(v) = max(moveError(v),abs(plan(1)-du(k)));
        end
    end
end
fprintf('Najwieksze naruszenie ograniczen: %.3g; blad wskaznikow: %.3g\n',maxConstraintError,maxMetricError);
fprintf('Maks. blad ruchu QP (FGPC-SL): %s\n',mat2str(moveError.',3));
assert(maxConstraintError<1e-8 && maxMetricError<1e-8 && all(moveError<1e-6));

%% 4. Niezalezny bilans objetosci dla zapisanych sterowan
maxPlantError = 0;
for scenario = 1:size(saved,1)
    for v = 1:numel(lambdaSL)
        control = saved(scenario,v).u;
        dist = disturbances{scenario};
        volume = [P.A1*P.h0(1);P.C2*P.h0(2)^3];
        checkY = zeros(numel(control)+1,1); checkY(1) = P.h0(2);
        for k = 1:numel(control)
            if k>P.delay, delayedU = control(k-P.delay); else, delayedU = P.u0; end
            dt = P.Ts/5;
            for sub = 1:5
                a = volume_rhs(P,volume,delayedU,dist(k));
                b = volume_rhs(P,volume+dt*a/2,delayedU,dist(k));
                c = volume_rhs(P,volume+dt*b/2,delayedU,dist(k));
                d = volume_rhs(P,volume+dt*c,delayedU,dist(k));
                volume = volume+dt*(a+2*b+2*c+d)/6;
            end
            checkY(k+1) = (volume(2)/P.C2)^(1/3);
        end
        maxPlantError = max(maxPlantError,max(abs(checkY-saved(scenario,v).y)));
    end
end
fprintf('Maks. roznica symulacji wysokosci i objetosci: %.3g cm\n',maxPlantError);
assert(maxPlantError<1e-4);
fprintf('AUDYT CZESCI 3 (FGPC-SL): WSZYSTKIE TESTY ZALICZONE.\n');

function val = get_u(u,idx,lastKnown,prev,u0)
    if idx<1, val = u0; elseif idx<=lastKnown, val = u(idx); else, val = prev; end
end

function dv = volume_rhs(P,v,u,d)
    h1 = max(v(1)/P.A1,0);
    h2 = max(v(2)/P.C2,0)^(1/3);
    dv = [u+d-P.a1*sqrt(h1); P.a1*sqrt(h1)-P.a2*sqrt(h2)];
end
