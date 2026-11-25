function plot_model_responses(obj,initial,tuned,steps,truth,t)
    selected = [1,6,7,12];
    fig = figure('Visible','off','Color','w','Position',[100 100 980 680]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    for k = 1:numel(selected)
        j = selected(k);
        ax = nexttile(tiles); ax.Color = 'w';
        stairs(ax,t,truth(:,j),'LineWidth',1); hold(ax,'on');
        stairs(ax,t,initial(:,j),'LineStyle','--','LineWidth',1);
        stairs(ax,t,tuned(:,j),'LineStyle','-.','LineWidth',1);
        grid(ax,'on'); xlabel(ax,'t [s]'); ylabel(ax,'h_2 [cm]');
        title(ax,sprintf('Skok F_{1in} %+g, F_D %+g cm^3/s',steps(j,1),steps(j,2)));
        if k==1
            lg = legend(ax,{'obiekt nieliniowy','TS 5 początkowy','TS 5 po strojeniu'},'NumColumns',3);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'porownanie_modeli_TS.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    close(fig);
end
