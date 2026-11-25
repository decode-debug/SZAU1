function pdc(obj,tControl,traces,controlCases)
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 840]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    for j = 1:numel(traces)
        y = traces(j).y;
        ref = traces(j).ref;
        ax = nexttile(tiles); hold(ax,'on'); grid(ax,'on'); ax.Color = 'w';
        stairs(ax,tControl,y(:,1),'LineWidth',1);
        stairs(ax,tControl,y(:,2),'LineStyle','--','LineWidth',1);
        stairs(ax,tControl,y(:,3),'LineWidth',1);
        stairs(ax,tControl,y(:,4),'LineStyle','-.','LineWidth',1);
        stairs(ax,tControl,ref,'LineStyle',':','LineWidth',1);
        xlabel(ax,'t [s]'); ylabel(ax,'y [cm]');
        xline(ax,400,':','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(ax,3800,':','Color',[.45 .45 .45],'HandleVisibility','off');
        xlim(ax,[tControl(1),tControl(end)]);
        title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,1),controlCases(j,2)));
        if j==1
            lg = legend(ax,{'DMC cz. 1 (bez ogr.)','DMC ograniczony', ...
                'PDC','PDC strojony','y_{zad}'},'NumColumns',3);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'regulacja_porownanie.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    close(fig);

    % Wykres sygnalu sterujacego u dla tych samych 4 scenariuszy.
    figU = figure('Visible','off','Color','w','Position',[100 100 1100 840]);
    theme(figU,'light');
    tilesU = tiledlayout(figU,2,2,'TileSpacing','compact','Padding','compact');
    for j = 1:numel(traces)
        u = traces(j).u;
        uPlot = [u; u(end,:)];
        ax = nexttile(tilesU); hold(ax,'on'); grid(ax,'on'); ax.Color = 'w';
        stairs(ax,tControl,uPlot(:,1),'LineWidth',1);
        stairs(ax,tControl,uPlot(:,2),'LineStyle','--','LineWidth',1);
        stairs(ax,tControl,uPlot(:,3),'LineWidth',1);
        stairs(ax,tControl,uPlot(:,4),'LineStyle','-.','LineWidth',1);
        xlabel(ax,'t [s]'); ylabel(ax,'u [cm^3/s]');
        xline(ax,400,':','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(ax,3800,':','Color',[.45 .45 .45],'HandleVisibility','off');
        xlim(ax,[tControl(1),tControl(end)]);
        title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,1),controlCases(j,2)));
        if j==1
            lg = legend(ax,{'DMC cz. 1 (bez ogr.)','DMC ograniczony', ...
                'PDC','PDC strojony'},'NumColumns',4);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(figU,fullfile(obj.OutputDir,'regulacja_sterowanie.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    close(figU);
end
