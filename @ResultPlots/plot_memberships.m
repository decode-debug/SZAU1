function plot_memberships(obj,x,weights,counts,tunedWeights)
    fig = figure('Visible','off','Color','w');
    tiledlayout(2,2);
    for m = 1:numel(counts)
        nexttile;
        plot(x,weights{m},'LineWidth',1.2);
        title(sprintf('%d reguly',counts(m)));
        xlabel('h_2 [cm]'); ylabel('stopien przynaleznosci'); grid on;
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'funkcje_przynaleznosci.png'),'Resolution',160);
    close(fig);

    fig = figure('Visible','off','Color','w');
    plot(x,tunedWeights,'LineWidth',1.5);
    title('Funkcje przynaleznosci po strojeniu modelu TS');
    xlabel('h_2 [cm]'); ylabel('stopien przynaleznosci'); grid on;
    exportgraphics(fig,fullfile(obj.OutputDir,'funkcje_przynaleznosci_strojone.png'), ...
                   'Resolution',160);
    close(fig);
end
