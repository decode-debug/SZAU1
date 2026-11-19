function plot_tuning(obj,t,traces,slTraces,cases,refStart,distStart)
    [~,j] = max(cases(:,1));
    fig = figure('Visible','off','Color','w','Position',[100 100 990 430]);
    tiles = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    masks = {t>=refStart+600 & t<=distStart-600, t>=distStart-100 & t<=distStart+2000};
    headings = {'Dojście do wartości zadanej (powiększenie)','Odpowiedź na skok zakłócenia'};
    for panel = 1:2
        ax = nexttile(tiles); ax.Toolbar.Visible = 'off';
        m = masks{panel};
        plot(ax,t(m),slTraces(j,2).y(m),'--','Color',[.35 .38 .41],'LineWidth',1.6);
        hold(ax,'on');
        plot(ax,t(m),slTraces(j,3).y(m),'Color',[0 .34 .64],'LineWidth',2);
        stairs(ax,t(m),traces(j).ref(m),':','Color',[.12 .47 .23],'LineWidth',1.6);
        if panel==2
            xline(ax,distStart,':','skok F_D','HandleVisibility','off');
        end
        xlabel(ax,'t [s]'); ylabel(ax,'h_2 [cm]'); title(ax,headings{panel});
        grid(ax,'on'); xlim(ax,[min(t(m)),max(t(m))]);
        if panel==1
            lg = legend(ax,{'SL przed strojeniem','SL po strojeniu','zadana'},'NumColumns',3);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'strojenie_SL.png'),'Resolution',160);
    close(fig);
end
