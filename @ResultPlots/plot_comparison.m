function plot_comparison(obj,t,traces,slTraces,cases,lambdaSL,selectedVariant,refStart,distStart)
    figY = figure('Visible','off','Color','w','Position',[100 100 1100 840]);
    theme(figY,'light');
    ty = tiledlayout(figY,2,2,'TileSpacing','compact','Padding','compact');
    figU = figure('Visible','off','Color','w','Position',[100 100 1100 840]);
    theme(figU,'light');
    tu = tiledlayout(figU,2,2,'TileSpacing','compact','Padding','compact');

    % Porownanie: PDC, PDC strojony, FGPC-SL (lam=1), FGPC-SL (wybrany lam=0.5), FGPC-SL (lam=0.2).
    slIdx = [2, selectedVariant, 4];
    names = ["PDC"; "PDC strojony"; ...
        sprintf("FGPC-SL (\\lambda = %g)",lambdaSL(slIdx(1))); ...
        sprintf("FGPC-SL (\\lambda = %g)",lambdaSL(slIdx(2))); ...
        sprintf("FGPC-SL (\\lambda = %g)",lambdaSL(slIdx(3)))];
    styles = {'--','-.',':','-','--'};

    for j = 1:numel(traces)
        outputs = [traces(j).y(:,3:4), slTraces(j,slIdx(1)).y, slTraces(j,slIdx(2)).y, slTraces(j,slIdx(3)).y];
        controls = [traces(j).u(:,3:4), slTraces(j,slIdx(1)).u, slTraces(j,slIdx(2)).u, slTraces(j,slIdx(3)).u];

        ay = nexttile(ty); ay.Color = 'w';
        linesY = stairs(ay,t,outputs); hold(ay,'on');
        refLine = stairs(ay,t,traces(j).ref,'LineStyle',':','LineWidth',1);
        ylabel(ay,'y [cm]');

        au = nexttile(tu); au.Color = 'w';
        linesU = stairs(au,t,[controls;controls(end,:)]);
        ylabel(au,'u [cm^3/s]');

        for m = 1:numel(names)
            set([linesY(m),linesU(m)],'LineStyle',styles{m},'LineWidth',1);
        end
        for ax = [ay,au]
            ax.Toolbar.Visible = 'off';
            grid(ax,'on'); xlabel(ax,'t [s]'); xlim(ax,[t(1),t(end)]);
            title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',cases(j,1),cases(j,2)));
            xline(ax,refStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
            xline(ax,distStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        end
        if j==1
            ly = legend(ay,[linesY(:);refLine],[names;"y_{zad}"],'NumColumns',3);
            lu = legend(au,linesU,names,'NumColumns',3);
            ly.Layout.Tile = 'south'; lu.Layout.Tile = 'south';
        end
    end
    exportgraphics(figY,fullfile(obj.OutputDir,'porownanie_wyjsc.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    exportgraphics(figU,fullfile(obj.OutputDir,'porownanie_sterowan.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    close(figY); close(figU);
end
