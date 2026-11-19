function steps(obj,t,responses,amplitudes,signal,filename)
    fig = figure('Visible','off','Color','w'); hold on;
    colors = lines(numel(amplitudes));
    for j = 1:numel(amplitudes)
        plot(t,responses(:,2*j-1),'-','Color',colors(j,:),'LineWidth',1.4,'DisplayName',sprintf('%+g nieliniowy',amplitudes(j)));
        plot(t,responses(:,2*j),'--','Color',colors(j,:),'LineWidth',1.2,'DisplayName',sprintf('%+g liniowy',amplitudes(j)));
    end
    xlabel('t [s]'); ylabel('h_2 [cm]'); title(['Skoki ' signal ' z punktu pracy']);
    grid on; legend('Location','eastoutside');
    exportgraphics(fig,fullfile(obj.OutputDir,filename),'Resolution',160); close(fig);
end
