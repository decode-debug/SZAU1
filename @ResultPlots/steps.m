function steps(obj,t,responses,amplitudes,signal,filename)
    fig = figure('Visible','off','Color','w'); theme(fig,'light'); hold on;
    for j = 1:numel(amplitudes)
        stairs(t,responses(:,2*j-1),'-','LineWidth',1,'DisplayName',sprintf('%+g nieliniowy',amplitudes(j)));
        stairs(t,responses(:,2*j),'--','LineWidth',1,'DisplayName',sprintf('%+g liniowy',amplitudes(j)));
    end
    xlabel('t [s]'); ylabel('y [cm]'); title(['Skoki ' signal ' z punktu pracy']);
    grid on; legend('Location','eastoutside');
    exportgraphics(fig,fullfile(obj.OutputDir,filename),'Resolution',160,'BackgroundColor','white'); close(fig);
end
