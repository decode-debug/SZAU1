function dmc(obj,tc,traces,scenarios)
    outdir = obj.OutputDir;
    fig = figure('Visible','off','Color','w','Position',[100 100 1050 1000]);
    tiledlayout(4,2);
    for j = 1:numel(traces)
        yNL = traces(j).yNL; yLin = traces(j).yLin;
        uNL = traces(j).uNL; uLin = traces(j).uLin;
        ref = traces(j).ref; dist = traces(j).dist;
        deltaRef = scenarios(j,1); deltaDist = scenarios(j,2);
        outputTile = 4*floor((j-1)/2)+mod(j-1,2)+1;
        nexttile(outputTile);
        plot(tc,yNL,'b',tc,yLin,'r--',tc,ref,'k:','LineWidth',1.2);
        grid on; ylabel('y(k) [cm]');
        title(sprintf('Skok y_{zad}: %+g cm, F_D: %+g cm^3/s',deltaRef,deltaDist));
        if j==1
            legend('obiekt nieliniowy','obiekt liniowy','y_{zad}(k)', ...
                   'Location','best');
        end

        nexttile(outputTile+2);
        stairs(tc,[uNL;uNL(end)],'b','LineWidth',1.2);
        hold on;
        stairs(tc,[uLin;uLin(end)],'r--','LineWidth',1.2);
        grid on; xlabel('t [s]'); ylabel('u(k) [cm^3/s]');
        plot_dmc_scenario(tc,yNL,yLin,uNL,uLin,ref,deltaRef,deltaDist,j,outdir);

        if j==2
            % Zachowaj dotychczasowy wykres dla pierwotnego testu.
            plot_dmc_original(tc,yNL,uNL,ref,dist,outdir);
        end
    end
    exportgraphics(fig,fullfile(outdir,'regulacja_DMC_porownanie.png'), ...
                   'Resolution',160);
    close(fig);
end

function plot_dmc_scenario(t,yNL,yLin,uNL,uLin,ref,deltaRef,deltaDist,j,outdir)
    % Wyjscie i sterowanie sa nierozlaczna para dla kazdej proby.
    fig = figure('Visible','off','Color','w','Position',[100 100 1000 650]);
    tiledlayout(2,1);
    ay = nexttile;
    plot(t,yNL,'b',t,yLin,'r--','LineWidth',1.2);
    hold on; stairs(t,ref,'k:','LineWidth',1.2);
    grid on; ylabel('y(k) [cm]');
    title(sprintf('Skok y_{zad}: %+g cm, F_D: %+g cm^3/s',deltaRef,deltaDist));
    legend('obiekt nieliniowy','obiekt liniowy','y_{zad}(k)','Location','best');
    au = nexttile;
    stairs(t,[uNL;uNL(end)],'b','LineWidth',1.2);
    hold on; stairs(t,[uLin;uLin(end)],'r--','LineWidth',1.2);
    grid on; xlabel('t [s]'); ylabel('u(k) [cm^3/s]');
    linkaxes([ay,au],'x');
    stem = sprintf('DMC_scenariusz_%d',j);
    exportgraphics(fig,fullfile(outdir,[stem '.png']),'Resolution',160);
    exportgraphics(fig,fullfile(outdir,[stem '.pdf']),'ContentType','vector');
    close(fig);
end

function plot_dmc_original(t,y,u,ref,dist,outdir)
    fig = figure('Visible','off','Color','w');
    tiledlayout(2,1);
    nexttile;
    plot(t,y,'b',t,ref,'k--','LineWidth',1.4);
    grid on; ylabel('y(k) [cm]');
    legend('obiekt nieliniowy','y_{zad}(k)','Location','best');
    title('DMC: skok y_{zad} +4 cm, skok F_D +20 cm^3/s');
    nexttile;
    stairs(t(1:end-1),u,'b','LineWidth',1.3);
    hold on; yyaxis right;
    stairs(t,dist,'r--','LineWidth',1.1);
    grid on; xlabel('t [s]'); ylabel('F_D [cm^3/s]');
    yyaxis left; ylabel('u(k) [cm^3/s]');
    ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
    exportgraphics(fig,fullfile(outdir,'regulacja_DMC.png'),'Resolution',160);
    close(fig);
end
