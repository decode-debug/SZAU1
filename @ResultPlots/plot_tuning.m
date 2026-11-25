function plot_tuning(obj,t,traces,slTraces,cases,lambdaSL,refStart,distStart)
    % Wplyw parametru lambda w FGPC-SL dla duzego skoku dodatniego (+8 cm) i ujemnego (-8 cm).
    [~,jPos] = max(cases(:,1));
    [~,jNeg] = min(cases(:,1));
    fig = figure('Visible','off','Color','w','Position',[100 100 1060 720]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    % Wykres dla wybranych wartosci lambda (2, 0.5, 0.2).
    showIdx = [1, 3, 4];
    styles = {'--','-','-.'};
    labels = arrayfun(@(l) sprintf('FGPC-SL (\\lambda = %g)',l),lambdaSL(showIdx),'UniformOutput',false);

    picked = [jPos, jNeg];
    for col = 1:2
        j = picked(col);
        % Gorny wiersz: wyjscie y
        ay = nexttile(tiles,col); hold(ay,'on'); ay.Color = 'w';
        for idx = 1:numel(showIdx)
            v = showIdx(idx);
            stairs(ay,t,slTraces(j,v).y,'LineStyle',styles{idx},'LineWidth',1);
        end
        stairs(ay,t,traces(j).ref,'LineStyle',':','LineWidth',1);
        xline(ay,refStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(ay,distStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        grid(ay,'on'); xlim(ay,[t(1),t(end)]);
        ylabel(ay,'y [cm]'); xlabel(ay,'t [s]');
        title(ay,sprintf('Wyjście y: zadanie %+g cm, F_D %+g cm^3/s',cases(j,1),cases(j,2)));

        % Dolny wiersz: sterowanie u
        au = nexttile(tiles,col+2); hold(au,'on'); au.Color = 'w';
        for idx = 1:numel(showIdx)
            v = showIdx(idx);
            u = slTraces(j,v).u;
            stairs(au,t,[u;u(end)],'LineStyle',styles{idx},'LineWidth',1);
        end
        xline(au,refStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(au,distStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
        grid(au,'on'); xlim(au,[t(1),t(end)]);
        ylabel(au,'u [cm^3/s]'); xlabel(au,'t [s]');
        title(au,sprintf('Sterowanie u: zadanie %+g cm, F_D %+g cm^3/s',cases(j,1),cases(j,2)));

        if col==1
            lg = legend(ay,[labels,{'y_{zad}'}],'NumColumns',4);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'strojenie_SL.pdf'), ...
        'ContentType','vector','BackgroundColor','white');
    close(fig);
end
