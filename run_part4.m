%% Zadanie 4: FDMC-NPL na tym samym obiekcie i scenariuszach co FDMC-SL
clear; close all; clc;
root = fileparts(mfilename('fullpath'));
source = fullfile(root,'wyniki_czesc3','dane_czesc3.mat');
if ~isfile(source), error('Uruchom najpierw run_part2 i run_part3.'); end
load(source,'P','ctrl','modelTuned','localSteps','traces','slTraces','controlCases', ...
    'tControl','refStart','distStart','predictionHorizons','horizons','lambdaSL','referenceTime');
outdir = fullfile(root,'wyniki_czesc4');
if ~isfolder(outdir), mkdir(outdir); end
options = optimoptions('quadprog','Algorithm','active-set','Display','off');
% Porownanie przy tych samych nastawach co dostrojony SL.
v = numel(horizons); ctrl.N = predictionHorizons(v); ctrl.lambda = lambdaSL(v);
Nu = horizons(v); Tf = referenceTime(v);
names = ["DMC cz. 1";"DMC ograniczony";"PDC";"PDC strojony"; ...
    "SL o nastawach PDC";"SL przed strojeniem";"SL po strojeniu";"NPL"];
task = tControl>=refStart & tControl<distStart; loadStep = tControl>=distStart;
rows = zeros(numel(traces)*numel(names),9);
for j=1:numel(traces)
    ref = traces(j).ref;
    [y,u,dt] = simulate_fdmc(P,localSteps,modelTuned.centers,ctrl,Nu,Tf, ...
        ref,traces(j).dist,options,modelTuned);
    results(j).y = y; results(j).u = u; results(j).decisionTime = dt;
    Y = [traces(j).y,cat(2,slTraces(j,:).y),y];
    U = [traces(j).u,cat(2,slTraces(j,:).u),u];
    for m=1:numel(names)
        e = Y(:,m)-ref; du = diff([P.u0;U(:,m)]);
        z = sign(controlCases(j,1))*Y(task,m);
        rows((j-1)*numel(names)+m,:) = [j,m,P.Ts*sum(abs(e(task))), ...
            P.Ts*sum(abs(e(loadStep))),max(cummax(z)-z),max(abs(e(loadStep))), ...
            max(abs(du)),sum(abs(du)),max(U(:,m))];
    end
end
quality = array2table(rows,'VariableNames',{'Scenariusz','ID','IAE_zadanie','IAE_zaklocenie', ...
    'Cofniecie_cm','Maks_blad_po_FD_cm','Maks_zmiana_u','Suma_zmian_u','Maks_u'});
quality.Regulator = names(quality.ID);
totals = table(names,accumarray(quality.ID,quality.IAE_zadanie), ...
    accumarray(quality.ID,quality.IAE_zaklocenie),'VariableNames',{'Regulator','IAE_zadanie','IAE_zaklocenie'});
totals.IAE = totals.IAE_zadanie+totals.IAE_zaklocenie;
disp(totals);
writetable(quality,fullfile(outdir,'porownanie_NPL.csv'));
writetable(totals,fullfile(outdir,'podsumowanie_NPL.csv'));
save(fullfile(outdir,'dane_czesc4.mat'),'results','tControl','controlCases','ctrl','Nu','Tf','quality','totals');
for signal=["y","u"]
    fig=figure('Visible','off','Color','w','Position',[100 100 1050 760]);
    theme(fig,'light');
    tiles=tiledlayout(fig,2,2,'TileSpacing','compact');
    for j=1:numel(traces)
        ax=nexttile(tiles); hold(ax,'on'); grid(ax,'on');
        ax.Color='w';
        curves=[traces(j).(signal)(:,1),traces(j).(signal)(:,4),slTraces(j,v).(signal),results(j).(signal)];
        if signal=="y"
            plot(ax,tControl,curves,'LineWidth',1.3);
            stairs(ax,tControl,traces(j).ref,'k:'); ylabel(ax,'h_2 [cm]');
            labels=[names([1,4,7,8]);"zadana"];
        else
            stairs(ax,tControl,[curves;curves(end,:)],'LineWidth',1.3);
            ylabel(ax,'F_{1in} [cm^3/s]'); labels=names([1,4,7,8]);
        end
        xline(ax,distStart,':','skok F_D','HandleVisibility','off');
        title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,:)));
        xlabel(ax,'t [s]');
        if j==1, lg=legend(ax,labels,'NumColumns',numel(labels)); lg.Layout.Tile='south'; end
    end
    exportgraphics(fig,fullfile(outdir,"NPL_"+signal+".png"),'Resolution',160,'BackgroundColor','white'); close(fig);
end
