function plot_comparison(obj,t,traces,slTraces,cases,names,refStart,distStart)
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
    exportgraphics(figY,fullfile(obj.OutputDir,'porownanie_wyjsc.png'),'Resolution',160);
    exportgraphics(figU,fullfile(obj.OutputDir,'porownanie_sterowan.png'),'Resolution',160);
    close(figY); close(figU);
end
