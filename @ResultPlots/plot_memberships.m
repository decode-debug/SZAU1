function plot_memberships(obj,x,weights,counts,tunedWeights)
    fig = figure('Visible','off','Color','w','Position',[100 100 980 680]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    for m = 1:numel(counts)
        ax = nexttile(tiles); ax.Color = 'w';
        plot(ax,x,weights{m},'LineWidth',1);
        if counts(m)==5, title(ax,'5 reguł'); else, title(ax,sprintf('%d reguły',counts(m))); end
        xlabel(ax,'y [cm]'); ylabel(ax,'stopień przynależności \mu_i(y)'); grid(ax,'on');
        ylim(ax,[0 1.05]);
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'funkcje_przynaleznosci.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    close(fig);

    fig = figure('Visible','off','Color','w','Position',[100 100 760 420]);
    theme(fig,'light');
    ax = axes(fig); ax.Color = 'w';
    plot(ax,x,tunedWeights,'LineWidth',1);
    title(ax,'Funkcje przynależności po strojeniu modelu TS (5 reguł)');
    xlabel(ax,'y [cm]'); ylabel(ax,'stopień przynależności \mu_i(y)'); grid(ax,'on');
    ylim(ax,[0 1.05]);
    exportgraphics(fig,fullfile(obj.OutputDir,'funkcje_przynaleznosci_strojone.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    close(fig);
end
