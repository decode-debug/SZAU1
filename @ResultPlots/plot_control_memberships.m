function plot_control_memberships(obj,x,initialWeights,tunedWeights,gamma)
    fig = figure('Visible','off','Color','w');
    tiledlayout(1,2);
    for panel = 1:2
        nexttile; hold on;
        if panel==1
            values = initialWeights;
        else
            values = tunedWeights;
        end
        plot(x,values,'LineWidth',1.3);
        if panel==1, title('PDC: funkcje jak w modelu');
        else, title(sprintf('Po strojeniu: gamma = %.2f',gamma)); end
        xlabel('h_2 [cm]'); ylabel('waga reguly'); grid on;
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'strojenie_regulatora.png'),'Resolution',160);
    close(fig);
end
