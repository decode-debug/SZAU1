function plot_model_responses(obj,initial,tuned,steps,truth,t)
    % Cztery reprezentatywne odpowiedzi; pelne RMSE liczone jest dla 12 prob.
    selected = [1,6,7,12];
    fig = figure('Visible','off','Color','w');
    tiledlayout(2,2);
    for k = 1:numel(selected)
        j = selected(k);
        yInitial = initial(:,j);
        yTuned = tuned(:,j);
        nexttile;
        plot(t,truth(:,j),'k',t,yInitial,'b--',t,yTuned,'r-.', ...
             'LineWidth',1.2);
        grid on; xlabel('t [s]'); ylabel('h_2 [cm]');
        title(sprintf('Skok F1in %+g, FD %+g',steps(j,1),steps(j,2)));
        if k==1
            legend('nieliniowy','TS poczatkowy','TS po strojeniu', ...
                   'Location','best');
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'porownanie_modeli_TS.png'),'Resolution',160);
    close(fig);
end
