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
    colors = [0.48 0.50 0.53; 0.85 0.33 0.10; 0.49 0.18 0.56; 0 0.45 0.74; 0.93 0.69 0.13];
    styles = {'--','-.',':','-','--'};
    widths = [1.15, 1.35, 1.35, 1.85, 1.3];

    for j = 1:numel(traces)
        outputs = [traces(j).y(:,3:4), slTraces(j,slIdx(1)).y, slTraces(j,slIdx(2)).y, slTraces(j,slIdx(3)).y];
        controls = [traces(j).u(:,3:4), slTraces(j,slIdx(1)).u, slTraces(j,slIdx(2)).u, slTraces(j,slIdx(3)).u];

        ay = nexttile(ty); ay.Color = 'w';
        linesY = stairs(ay,t,outputs); hold(ay,'on');
        refLine = stairs(ay,t,traces(j).ref,':','Color',[0.12 0.47 0.23],'LineWidth',1.5);
        ylabel(ay,'h_2 [cm]');

        au = nexttile(tu); au.Color = 'w';
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
    exportgraphics(figY,fullfile(obj.OutputDir,'porownanie_wyjsc.png'), ...
        'Resolution',160,'BackgroundColor','white');
    exportgraphics(figU,fullfile(obj.OutputDir,'porownanie_sterowan.png'), ...
        'Resolution',160,'BackgroundColor','white');
    close(figY); close(figU);
end
