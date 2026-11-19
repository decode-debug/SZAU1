function plot_all_models(obj,linearOutputs,modelOutputs,counts,steps,truth,t)
    % Porownanie odpowiedzi 2, 3, 4 i 5 regul z liniowym punktem odniesienia.
    selected = [1,6,7,12];
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 850]);
    tiledlayout(2,2);
    colors = lines(numel(counts));
    for tile = 1:numel(selected)
        j = selected(tile);
        nexttile; hold on;
        plot(t,truth(:,j),'k','LineWidth',1.7,'DisplayName','nieliniowy');
        yLinear = linearOutputs(:,j);
        plot(t,yLinear,'k:','LineWidth',1.2,'DisplayName','liniowy');
        for m = 1:numel(counts)
            yTS = modelOutputs(:,j,m);
            plot(t,yTS,'Color',colors(m,:),'LineWidth',1.1, ...
                 'DisplayName',sprintf('TS %d',counts(m)));
        end
        grid on; xlabel('t [s]'); ylabel('h_2 [cm]');
        title(sprintf('Skok F1in %+g, FD %+g',steps(j,1),steps(j,2)));
        if tile==1, legend('Location','best'); end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'modele_2_5.png'),'Resolution',160);
    close(fig);
end
