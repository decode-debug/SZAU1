function pdc(obj,tControl,traces,controlCases)
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 840]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    for j = 1:numel(traces)
        y = traces(j).y;
        ref = traces(j).ref;
        ax = nexttile(tiles); hold(ax,'on'); grid(ax,'on'); ax.Color = 'w';
        stairs(ax,tControl,y(:,1),'Color',[0.48 0.48 0.48],'LineWidth',1.2);
        stairs(ax,tControl,y(:,2),'k--','LineWidth',1.2);
        stairs(ax,tControl,y(:,3),'Color',[0 0.45 0.74],'LineWidth',1.35);
        stairs(ax,tControl,y(:,4),'Color',[0.85 0.33 0.10],'LineStyle','-.','LineWidth',1.4);
        stairs(ax,tControl,ref,'Color',[0.12 0.47 0.23],'LineStyle',':','LineWidth',1.4);
        xlabel(ax,'t [s]'); ylabel(ax,'h_2 [cm]');
        xlim(ax,[tControl(1),tControl(end)]);
        title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,1),controlCases(j,2)));
        if j==1
            lg = legend(ax,{'DMC cz. 1 (bez ogr.)','Liniowy DMC (ogr.)', ...
                'Rozmyty PDC','Rozmyty PDC strojony','zadana'},'NumColumns',3);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'regulacja_porownanie.png'), ...
        'Resolution',160,'BackgroundColor','white');
    close(fig);
end
