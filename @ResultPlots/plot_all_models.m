function plot_all_models(obj,linearOutputs,modelOutputs,counts,steps,truth,t)
    selected = [1,6,7,12];
    fig = figure('Visible','off','Color','w','Position',[100 100 1080 780]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    colors = lines(numel(counts));
    for tile = 1:numel(selected)
        j = selected(tile);
        ax = nexttile(tiles); hold(ax,'on'); ax.Color = 'w';
        plot(ax,t,truth(:,j),'k','LineWidth',1.7,'DisplayName','obiekt nieliniowy');
        plot(ax,t,linearOutputs(:,j),'k:','LineWidth',1.3,'DisplayName','liniowy');
        for m = 1:numel(counts)
            plot(ax,t,modelOutputs(:,j,m),'Color',colors(m,:),'LineWidth',1.2, ...
                'DisplayName',sprintf('TS %d',counts(m)));
        end
        grid(ax,'on'); xlabel(ax,'t [s]'); ylabel(ax,'h_2 [cm]');
        title(ax,sprintf('Skok F_{1in} %+g, F_D %+g cm^3/s',steps(j,1),steps(j,2)));
        if tile==1
            lg = legend(ax,'NumColumns',6);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'modele_2_5.png'), ...
        'Resolution',160,'BackgroundColor','white');
    close(fig);
end
