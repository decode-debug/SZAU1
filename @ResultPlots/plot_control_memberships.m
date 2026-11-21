function plot_control_memberships(obj,x,initialWeights,tunedWeights,gamma)
    fig = figure('Visible','off','Color','w','Position',[100 100 960 390]);
    theme(fig,'light');
    tiles = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    for panel = 1:2
        ax = nexttile(tiles); hold(ax,'on'); ax.Color = 'w';
        if panel==1
            plot(ax,x,initialWeights,'LineWidth',1.35);
            title(ax,'PDC: funkcje przynależności jak w modelu');
        else
            plot(ax,x,tunedWeights,'LineWidth',1.35);
            title(ax,sprintf('Po strojeniu regulatora: \\gamma = %.2f',gamma));
        end
        xlabel(ax,'h_2 [cm]'); ylabel(ax,'znormalizowana waga reguły'); grid(ax,'on');
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'strojenie_regulatora.png'), ...
        'Resolution',160,'BackgroundColor','white');
    close(fig);
end
