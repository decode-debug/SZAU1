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
if ~isfolder(outdir), mkdir(outdir); end

set(groot,'defaultFigureColor','w','defaultAxesColor','w', ...
    'defaultAxesXColor','k','defaultAxesYColor','k', ...
    'defaultTextColor','k','defaultLegendColor','w', ...
    'defaultLegendTextColor','k');

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
        [y,u,decisionTime] = simulate_sl(P,localSteps,modelTuned.centers, ...
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

plot_comparison(tControl,traces,slTraces,controlCases,methodNames,refStart,distStart,outdir);
plot_tuning(tControl,traces,slTraces,controlCases,refStart,distStart,outdir);
writetable(slSummary,fullfile(outdir,'wskazniki_SL.csv'));
writetable(comparisonSummary,fullfile(outdir,'porownanie_regulatorow.csv'));
writetable(validationSummary,fullfile(outdir,'walidacja_SL.csv'));
save(fullfile(outdir,'dane_czesc3.mat'),'P','modelTuned','ctrl','controlCases', ...
    'tControl','traces','slTraces','horizons','predictionHorizons','lambdaSL', ...
    'referenceTime','refStart','distStart','localSteps','slSummary','comparisonSummary', ...
    'totalIAE','validationCases','validationTraces','validationSummary');

%% Funkcje pomocnicze
function [y,u,decisionTime] = simulate_sl(P,localSteps,centers,ctrl,Nu,Tf,reference,disturbance,options)
    n = numel(reference)-1;
    y = zeros(n+1,1); y(1) = P.h0(2);
    u = zeros(n,1);
    decisionTime = zeros(n,1);
    h = P.h0;
    pastDu = zeros(ctrl.D-1,1);
    cumulative = tril(ones(Nu));
    Aineq = [cumulative;-cumulative];
    lb = -ctrl.duMax*ones(Nu,1);
    ub = ctrl.duMax*ones(Nu,1);
    target = reference(1);
    alpha = 0;
    if Tf>0, alpha = exp(-P.Ts/Tf); end

    for k = 1:n
        if k==1, previousU=P.u0; else, previousU=u(k-1); end
        started = tic;

        % (4.97): kazda rzedna jest suma wazona lokalnych rzednych skoku.
        weights = ts_memberships(y(k),centers);
        s = localSteps*weights.';
        [M,freeMatrix] = dmc_matrices(s,ctrl.N,Nu);
        free = y(k)+freeMatrix*pastDu;

        H = 2*(M'*M+ctrl.lambda*eye(Nu));
        % Filtr przyczynowy. Regulator zna tylko biezaca wartosc zadana.
        target = alpha*target+(1-alpha)*reference(k);
        f = -2*M'*(target-free);
        bineq = [(ctrl.uMax-previousU)*ones(Nu,1); ...
                 (previousU-ctrl.uMin)*ones(Nu,1)];
        [moves,~,exitflag] = quadprog(H,f,Aineq,bineq,[],[], ...
                                      lb,ub,zeros(Nu,1),options);
        if exitflag<=0
            error('quadprog nie rozwiazal probki k=%d (exitflag=%d).',k,exitflag);
        end
        if any(Aineq*moves>bineq+1e-8) || any(moves<lb-1e-8 | moves>ub+1e-8)
            error('Plan sterowania narusza ograniczenia w probce k=%d.',k);
        end
        decisionTime(k) = toc(started);

        du = moves(1);
        u(k) = previousU+du;
        pastDu = [du;pastDu(1:end-1)];
        if k>P.delay, actualU=u(k-P.delay); else, actualU=P.u0; end
        h = rk4_step(@(state) tank_rhs(P,state,actualU,disturbance(k)), ...
                     h,P.Ts,5);
        h = max(h,0.05);
        y(k+1) = h(2);
    end
end

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

function dh = tank_rhs(P,h,u,d)
    f2 = P.a1*sqrt(max(h(1),0));
    f3 = P.a2*sqrt(max(h(2),0));
    cap2 = 3*P.C2*max(h(2),0.05)^2;
    dh = [(u+d-f2)/P.A1;(f2-f3)/cap2];
end

function state = rk4_step(rhs,state,T,n)
    dt = T/n;
    for j = 1:n
        a = rhs(state);
        b = rhs(state+dt*a/2);
        c = rhs(state+dt*b/2);
        d = rhs(state+dt*c);
        state = state+dt*(a+2*b+2*c+d)/6;
    end
end

function plot_comparison(t,traces,slTraces,cases,names,refStart,distStart,outdir)
    figY = figure('Visible','off','Color','w','Position',[100 100 1100 850]);
    ty = tiledlayout(figY,2,2,'TileSpacing','compact','Padding','compact');
    figU = figure('Visible','off','Color','w','Position',[100 100 1100 850]);
    tu = tiledlayout(figU,2,2,'TileSpacing','compact','Padding','compact');
    colors = [0.44 0.47 0.50; 0.79 0.33 0.05; 0.52 0.27 0.64; 0.18 0.23 0.27; 0 0.34 0.64];
    styles = {'--','-.',':','--','-'};
    widths = [1.1,1.3,1.4,1.4,1.8];
    for j = 1:numel(traces)
        outputs = [traces(j).y(:,3:4),cat(2,slTraces(j,:).y)];
        controls = [traces(j).u(:,3:4),cat(2,slTraces(j,:).u)];
        ay = nexttile(ty);
        linesY = plot(ay,t,outputs); hold(ay,'on');
        refLine = stairs(ay,t,traces(j).ref,':','Color',[0.12 0.47 0.23],'LineWidth',1.5);
        ylabel(ay,'h_2 [cm]');
        au = nexttile(tu);
        linesU = stairs(au,t,[controls;controls(end,:)]);
        ylabel(au,'F_{1in} [cm^3/s]');
        for m = 1:numel(names)
            set([linesY(m),linesU(m)],'LineStyle',styles{m},'Color',colors(m,:),'LineWidth',widths(m));
        end
        for ax = [ay,au]
            ax.Toolbar.Visible = 'off';
            grid(ax,'on'); xlabel(ax,'t [s]'); xlim(ax,[t(1),t(end)]);
            title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',cases(j,1),cases(j,2)));
            xline(ax,refStart,':','zadanie','Color',[.45 .45 .45],'HandleVisibility','off');
            xline(ax,distStart,':','skok F_D','Color',[.45 .45 .45],'HandleVisibility','off');
        end
        if j==1
            ly = legend(ay,[linesY(:);refLine],[names;"zadana"],'NumColumns',3);
            lu = legend(au,linesU,names,'NumColumns',3);
            ly.Layout.Tile = 'south'; lu.Layout.Tile = 'south';
        end
    end
    exportgraphics(figY,fullfile(outdir,'porownanie_wyjsc.png'),'Resolution',160);
    exportgraphics(figU,fullfile(outdir,'porownanie_sterowan.png'),'Resolution',160);
    close(figY); close(figU);
end

function plot_tuning(t,traces,slTraces,cases,refStart,distStart,outdir)
    [~,j] = max(cases(:,1));
    fig = figure('Visible','off','Color','w','Position',[100 100 990 430]);
    tiles = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    masks = {t>=refStart+600 & t<=distStart-600, t>=distStart-100 & t<=distStart+2000};
    headings = {'Dojście do wartości zadanej (powiększenie)','Odpowiedź na skok zakłócenia'};
    for panel = 1:2
        ax = nexttile(tiles); ax.Toolbar.Visible = 'off';
        m = masks{panel};
        plot(ax,t(m),slTraces(j,2).y(m),'--','Color',[.35 .38 .41],'LineWidth',1.6);
        hold(ax,'on');
        plot(ax,t(m),slTraces(j,3).y(m),'Color',[0 .34 .64],'LineWidth',2);
        stairs(ax,t(m),traces(j).ref(m),':','Color',[.12 .47 .23],'LineWidth',1.6);
        if panel==2
            xline(ax,distStart,':','skok F_D','HandleVisibility','off');
        end
        xlabel(ax,'t [s]'); ylabel(ax,'h_2 [cm]'); title(ax,headings{panel});
        grid(ax,'on'); xlim(ax,[min(t(m)),max(t(m))]);
        if panel==1
            lg = legend(ax,{'SL przed strojeniem','SL po strojeniu','zadana'},'NumColumns',3);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(outdir,'strojenie_SL.png'),'Resolution',160);
    close(fig);
end
