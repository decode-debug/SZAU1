function pdc(obj,tControl,traces,controlCases)
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 900]);
    tiledlayout(2,2);
    for j = 1:numel(traces)
        y = traces(j).y;
        ref = traces(j).ref;
        nexttile;
        plot(tControl,y(:,1),'Color',[0.5 0.5 0.5],'LineWidth',1.15);
        hold on;
        plot(tControl,y(:,2),'k--',tControl,y(:,3),'b', ...
             tControl,y(:,4),'r-.',tControl,ref,'g:', ...
             'LineWidth',1.15);
        grid on; xlabel('t [s]'); ylabel('h_2 [cm]');
        title(sprintf('Zadanie %+g cm, F_D %+g cm^3/s', ...
              controlCases(j,1),controlCases(j,2)));
        if j==1
            legend('DMC czesc 1','liniowy z ograniczeniami', ...
                'rozmyty PDC','rozmyty po strojeniu','zadana', ...
                'Location','best');
        end
    end
    exportgraphics(fig,fullfile(obj.OutputDir,'regulacja_porownanie.png'),'Resolution',160);
    close(fig);
end
