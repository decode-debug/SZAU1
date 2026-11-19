function npl(obj,tControl,traces,slTraces,v,results,names,controlCases,distStart)
    for signal=["y","u"]
        fig=figure('Visible','off','Color','w','Position',[100 100 1050 760]);
        theme(fig,'light');
        tiles=tiledlayout(fig,2,2,'TileSpacing','compact');
        for j=1:numel(traces)
            ax=nexttile(tiles); hold(ax,'on'); grid(ax,'on');
            ax.Color='w';
            curves=[traces(j).(signal)(:,1),traces(j).(signal)(:,4),slTraces(j,v).(signal),results(j).(signal)];
            if signal=="y"
                plot(ax,tControl,curves,'LineWidth',1.3);
                stairs(ax,tControl,traces(j).ref,'k:'); ylabel(ax,'h_2 [cm]');
                labels=[names([1,4,7,8]);"zadana"];
            else
                stairs(ax,tControl,[curves;curves(end,:)],'LineWidth',1.3);
                ylabel(ax,'F_{1in} [cm^3/s]'); labels=names([1,4,7,8]);
            end
            xline(ax,distStart,':','skok F_D','HandleVisibility','off');
            title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,:)));
            xlabel(ax,'t [s]');
            if j==1, lg=legend(ax,labels,'NumColumns',numel(labels)); lg.Layout.Tile='south'; end
        end
        exportgraphics(fig,fullfile(obj.OutputDir,"NPL_"+signal+".png"),'Resolution',160,'BackgroundColor','white'); close(fig);
    end
end
