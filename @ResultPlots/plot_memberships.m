function plot_memberships(obj,x,weights,counts,tunedWeights)
    fig = figure('Visible','off','Color','w','Position',[100 100 980 680]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    for m = 1:numel(counts)
        ax = nexttile(tiles); ax.Color = 'w';
        plot(ax,x,weights{m},'LineWidth',1.35);
        if counts(m)==5, title(ax,'5 reguł'); else, title(ax,sprintf('%d reguły',counts(m))); end
        xlabel(ax,'h_2 [cm]'); ylabel(ax,'stopień przynależności'); grid(ax,'on');
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'funkcje_przynaleznosci.png'), ...
        'Resolution',160,'BackgroundColor','white');
    close(fig);

    fig = figure('Visible','off','Color','w','Position',[100 100 760 420]);
    theme(fig,'light');
    ax = axes(fig); ax.Color = 'w';
    plot(ax,x,tunedWeights,'LineWidth',1.5);
    title(ax,'Funkcje przynależności po strojeniu modelu TS (5 reguł)');
    xlabel(ax,'h_2 [cm]'); ylabel(ax,'stopień przynależności'); grid(ax,'on');
    exportgraphics(fig,fullfile(obj.OutputDir,'funkcje_przynaleznosci_strojone.png'), ...
        'Resolution',160,'BackgroundColor','white');
    close(fig);
end
