function npl(obj,tControl,traces,slTraces,nplTraces,selectedVariant,lambdaNPL,controlCases,refStart,distStart)
    % 1-2. Glowne porownanie wszystkich pokolen regulatorow: DMC cz. 1, PDC strojony, FGPC-SL i FGPC-NPL.
    lamSel = lambdaNPL(selectedVariant);
    baseLabels = [ ...
        "DMC cz. 1"; ...
        "PDC strojony"; ...
        sprintf("FGPC-SL (\\lambda = %g)",lamSel); ...
        sprintf("FGPC-NPL (\\lambda = %g)",lamSel)];
    colors = [0.48 0.50 0.53; 0.85 0.33 0.10; 0.93 0.69 0.13; 0 0.45 0.74];
    styles = {'--','-.','-','-'};
    widths = [1.15, 1.35, 1.45, 1.85];

    for signal = ["y","u"]
        fig = figure('Visible','off','Color','w','Position',[100 100 1080 800]);
        theme(fig,'light');
        tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
        for j = 1:numel(traces)
            ax = nexttile(tiles); hold(ax,'on'); grid(ax,'on'); ax.Color = 'w';
            curves = [traces(j).(signal)(:,1), traces(j).(signal)(:,4), ...
                      slTraces(j,selectedVariant).(signal), nplTraces(j,selectedVariant).(signal)];
            if signal=="y"
                hL = plot(ax,tControl,curves);
                refL = stairs(ax,tControl,traces(j).ref,'Color',[0.12 0.47 0.23],'LineStyle',':','LineWidth',1.45);
                ylabel(ax,'h_2 [cm]');
            else
                hL = stairs(ax,tControl,[curves;curves(end,:)]);
                ylabel(ax,'F_{1in} [cm^3/s]');
            end
            for m = 1:numel(baseLabels)
                set(hL(m),'Color',colors(m,:),'LineStyle',styles{m},'LineWidth',widths(m));
            end
            xline(ax,refStart,':','zadanie','Color',[.45 .45 .45],'HandleVisibility','off');
            xline(ax,distStart,':','skok F_D','Color',[.45 .45 .45],'HandleVisibility','off');
            xlim(ax,[tControl(1),tControl(end)]);
            title(ax,sprintf('Zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,:)));
            xlabel(ax,'t [s]');
            if j==1
                if signal=="y"
                    lg = legend(ax,[hL(:);refL],[baseLabels;"zadana"],'NumColumns',5);
                else
                    lg = legend(ax,hL,baseLabels,'NumColumns',4);
                end
                lg.Layout.Tile = 'south';
            end
        end
        exportgraphics(fig,fullfile(obj.OutputDir,"NPL_"+signal+".png"), ...
            'Resolution',160,'BackgroundColor','white');
        close(fig);
    end

    % 3. Wplyw parametru lambda w FGPC-NPL (analogia do Rys. 4.45 ze skryptu) + porownanie z FGPC-SL.
    [~,jPos] = max(controlCases(:,1));
    [~,jNeg] = min(controlCases(:,1));
    fig = figure('Visible','off','Color','w','Position',[100 100 1060 720]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    lColors = [0.49 0.18 0.56; 0.47 0.67 0.19; 0 0.45 0.74; 0.85 0.33 0.10];
    lStyles = {'--','-.','-','-'};
    lWidths = [1.25, 1.35, 1.85, 1.35];
    lLabels = arrayfun(@(l) sprintf('FGPC-NPL (\\lambda = %g)',l),lambdaNPL,'UniformOutput',false);

    picked = [jPos, jNeg];
    for col = 1:2
        j = picked(col);
        ay = nexttile(tiles,col); hold(ay,'on'); ay.Color = 'w';
        for v = 1:numel(lambdaNPL)
            plot(ay,tControl,nplTraces(j,v).y,'Color',lColors(v,:),'LineStyle',lStyles{v},'LineWidth',lWidths(v));
        end
        stairs(ay,tControl,traces(j).ref,'k:','LineWidth',1.4);
        xline(ay,refStart,':','zadanie','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(ay,distStart,':','skok F_D','Color',[.45 .45 .45],'HandleVisibility','off');
        grid(ay,'on'); xlim(ay,[tControl(1),tControl(end)]);
        ylabel(ay,'h_2 [cm]'); xlabel(ay,'t [s]');
        title(ay,sprintf('Wyjście h_2 (NPL): zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,:)));

        au = nexttile(tiles,col+2); hold(au,'on'); au.Color = 'w';
        for v = 1:numel(lambdaNPL)
            u = nplTraces(j,v).u;
            stairs(au,tControl,[u;u(end)],'Color',lColors(v,:),'LineStyle',lStyles{v},'LineWidth',lWidths(v));
        end
        xline(au,refStart,':','zadanie','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(au,distStart,':','skok F_D','Color',[.45 .45 .45],'HandleVisibility','off');
        grid(au,'on'); xlim(au,[tControl(1),tControl(end)]);
        ylabel(au,'F_{1in} [cm^3/s]'); xlabel(au,'t [s]');
        title(au,sprintf('Sterowanie F_{1in} (NPL): zadanie %+g cm, F_D %+g cm^3/s',controlCases(j,:)));

        if col==1
            lg = legend(ay,[lLabels,{'zadana'}],'NumColumns',5);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'strojenie_NPL.png'), ...
        'Resolution',160,'BackgroundColor','white');
    close(fig);
end
