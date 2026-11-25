function npl(obj,tControl,traces,slTraces,nplTraces,selectedVariant,lambdaNPL,controlCases,refStart,distStart)
    % Porownanie regulatorow: DMC cz. 1, PDC strojony, FGPC-SL i FGPC-NPL.
    lamSel = lambdaNPL(selectedVariant);
    baseLabels = [ ...
        "DMC cz. 1"; ...
        "PDC strojony"; ...
        sprintf("FGPC-SL (\\lambda = %g)",lamSel); ...
        sprintf("FGPC-NPL (\\lambda = %g)",lamSel)];
    styles = {'--','-.','-','-'};

    for signal = ["y","u"]
        fig = figure('Visible','off','Color','w','Position',[100 100 1080 800]);
        theme(fig,'light');
        tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
        for j = 1:numel(traces)
            ax = nexttile(tiles); hold(ax,'on'); grid(ax,'on'); ax.Color = 'w';
            curves = [traces(j).(signal)(:,1), traces(j).(signal)(:,4), ...
                      slTraces(j,selectedVariant).(signal), nplTraces(j,selectedVariant).(signal)];
            if signal=="y"
                hL = stairs(ax,tControl,curves);
                refL = stairs(ax,tControl,traces(j).ref,'LineStyle',':','LineWidth',1);
                ylabel(ax,'y [cm]');
            else
                hL = stairs(ax,tControl,[curves;curves(end,:)]);
                ylabel(ax,'u [cm^3/s]');
            end
            for m = 1:numel(baseLabels)
                set(hL(m),'LineStyle',styles{m},'LineWidth',1);
            end
            xline(ax,refStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
            xline(ax,distStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
            xlim(ax,[tControl(1),tControl(end)]);
            title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,:)));
            xlabel(ax,'t [s]');
            if j==1
                if signal=="y"
                    lg = legend(ax,[hL(:);refL],[baseLabels;"y_{zad}"],'NumColumns',5);
                else
                    lg = legend(ax,hL,baseLabels,'NumColumns',4);
                end
                lg.Layout.Tile = 'south';
            end
        end
        exportgraphics(fig,fullfile(obj.OutputDir,"NPL_"+signal+".pdf"), ...
            'ContentType','vector','BackgroundColor','white');
        close(fig);
    end

    % Wplyw parametru lambda w FGPC-NPL.
    [~,jPos] = max(controlCases(:,1));
    [~,jNeg] = min(controlCases(:,1));
    fig = figure('Visible','off','Color','w','Position',[100 100 1060 720]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    % Wykres dla wybranych wartosci lambda (2, 0.5, 0.2).
    showIdx = [1, 3, 4];
    lStyles = {'--','-','-.'};
    lLabels = arrayfun(@(l) sprintf('FGPC-NPL (\\lambda = %g)',l),lambdaNPL(showIdx),'UniformOutput',false);

    picked = [jPos, jNeg];
    for col = 1:2
        j = picked(col);
        ay = nexttile(tiles,col); hold(ay,'on'); ay.Color = 'w';
        for idx = 1:numel(showIdx)
            v = showIdx(idx);
            stairs(ay,tControl,nplTraces(j,v).y,'LineStyle',lStyles{idx},'LineWidth',1);
        end
        stairs(ay,tControl,traces(j).ref,'LineStyle',':','LineWidth',1);
        xline(ay,refStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(ay,distStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        grid(ay,'on'); xlim(ay,[tControl(1),tControl(end)]);
        ylabel(ay,'y [cm]'); xlabel(ay,'t [s]');
        title(ay,sprintf('Wyjście y (NPL): zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,:)));

        au = nexttile(tiles,col+2); hold(au,'on'); au.Color = 'w';
        for idx = 1:numel(showIdx)
            v = showIdx(idx);
            u = nplTraces(j,v).u;
            stairs(au,tControl,[u;u(end)],'LineStyle',lStyles{idx},'LineWidth',1);
        end
        xline(au,refStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(au,distStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        grid(au,'on'); xlim(au,[tControl(1),tControl(end)]);
        ylabel(au,'u [cm^3/s]'); xlabel(au,'t [s]');
        title(au,sprintf('Sterowanie u (NPL): zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,:)));

        if col==1
            lg = legend(ay,[lLabels,{'y_{zad}'}],'NumColumns',4);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'strojenie_NPL.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    close(fig);
end
