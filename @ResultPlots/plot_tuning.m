function plot_tuning(obj,t,traces,slTraces,cases,lambdaSL,refStart,distStart)
    % Wplyw parametru lambda w FGPC-SL dla duzego skoku dodatniego (+8 cm) i ujemnego (-8 cm).
    [~,jPos] = max(cases(:,1));
    [~,jNeg] = min(cases(:,1));
    fig = figure('Visible','off','Color','w','Position',[100 100 1060 720]);
    theme(fig,'light');
    tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    colors = [0.49 0.18 0.56; 0.47 0.67 0.19; 0 0.45 0.74; 0.85 0.33 0.10];
    styles = {'--','-.','-','-'};
    widths = [1.25, 1.35, 1.85, 1.35];
    labels = arrayfun(@(l) sprintf('FGPC-SL (\\lambda = %g)',l),lambdaSL,'UniformOutput',false);

    picked = [jPos, jNeg];
    for col = 1:2
        j = picked(col);
        % Gorny wiersz: wyjscie h_2
        ay = nexttile(tiles,col); hold(ay,'on'); ay.Color = 'w';
        for v = 1:numel(lambdaSL)
            plot(ay,t,slTraces(j,v).y,'Color',colors(v,:),'LineStyle',styles{v},'LineWidth',widths(v));
        end
        stairs(ay,t,traces(j).ref,'k:','LineWidth',1.4);
        xline(ay,refStart,':','zadanie','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(ay,distStart,':','skok F_D','Color',[.45 .45 .45],'HandleVisibility','off');
        grid(ay,'on'); xlim(ay,[t(1),t(end)]);
        ylabel(ay,'h_2 [cm]'); xlabel(ay,'t [s]');
        title(ay,sprintf('Wyjście h_2: zadanie %+g cm, F_D %+g cm^3/s',cases(j,1),cases(j,2)));

        % Dolny wiersz: sterowanie F_1in
        au = nexttile(tiles,col+2); hold(au,'on'); au.Color = 'w';
        for v = 1:numel(lambdaSL)
            u = slTraces(j,v).u;
            stairs(au,t,[u;u(end)],'Color',colors(v,:),'LineStyle',styles{v},'LineWidth',widths(v));
        end
        xline(au,refStart,':','zadanie','Color',[.45 .45 .45],'HandleVisibility','off');
        xline(au,distStart,':','skok F_D','Color',[.45 .45 .45],'HandleVisibility','off');
        grid(au,'on'); xlim(au,[t(1),t(end)]);
        ylabel(au,'F_{1in} [cm^3/s]'); xlabel(au,'t [s]');
        title(au,sprintf('Sterowanie F_{1in}: zadanie %+g cm, F_D %+g cm^3/s',cases(j,1),cases(j,2)));

        if col==1
            lg = legend(ay,[labels,{'zadana'}],'NumColumns',5);
            lg.Layout.Tile = 'south';
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'strojenie_SL.png'), ...
        'Resolution',160,'BackgroundColor','white');
    close(fig);
end
