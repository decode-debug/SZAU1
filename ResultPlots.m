classdef ResultPlots
% Zbior metod generujacych wektorowe wykresy PDF dla wszystkich 4 czesci projektu.
    properties (SetAccess = private)
        OutputDir
    end

    methods
        function obj = ResultPlots(outdir)
            obj.OutputDir = outdir;
            if ~isfolder(outdir), mkdir(outdir); end
            set(groot,'defaultFigureColor','w','defaultAxesColor','w', ...
                'defaultAxesXColor','k','defaultAxesYColor','k', ...
                'defaultTextColor','k','defaultLegendColor','w', ...
                'defaultLegendTextColor','k');
        end

        function steps(obj,t,responses,amplitudes,signal,filename)
            fig = obj.new_figure([100 100 760 440]);
            ax = axes(fig); hold(ax,'on'); ax.Color = 'w';
            for j = 1:numel(amplitudes)
                stairs(ax,t,responses(:,2*j-1),'-','LineWidth',1, ...
                    'DisplayName',sprintf('%+g nieliniowy',amplitudes(j)));
                stairs(ax,t,responses(:,2*j),'--','LineWidth',1, ...
                    'DisplayName',sprintf('%+g liniowy',amplitudes(j)));
            end
            xlabel(ax,'t [s]'); ylabel(ax,'y [cm]');
            title(ax,['Skoki ' signal ' z punktu pracy']);
            grid(ax,'on'); legend(ax,'Location','eastoutside');
            obj.save_pdf(fig,filename);
        end

        function dmc(obj,tc,traces,scenarios,refStart,distStart)
            fig = obj.new_figure([100 100 1080 960]);
            tiles = tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');
            for j = 1:numel(traces)
                tr = traces(j);
                row = 4*floor((j-1)/2) + mod(j-1,2) + 1;
                ay = nexttile(tiles,row); hold(ay,'on'); ay.Color = 'w';
                stairs(ay,tc,tr.yNL,'-','LineWidth',1);
                stairs(ay,tc,tr.yLin,'--','LineWidth',1);
                stairs(ay,tc,tr.ref,':','LineWidth',1);
                obj.decorate_time_axis(ay,tc,refStart,distStart,'y [cm]', ...
                    sprintf('Zadanie %+g cm, F_D %+g cm^3/s',scenarios(j,1),scenarios(j,2)));
                if j==1
                    lg = legend(ay,{'obiekt nieliniowy','obiekt liniowy','y_{zad}'},'NumColumns',3);
                    lg.Layout.Tile = 'south';
                end

                au = nexttile(tiles,row+2); hold(au,'on'); au.Color = 'w';
                stairs(au,tc,[tr.uNL;tr.uNL(end)],'-','LineWidth',1);
                stairs(au,tc,[tr.uLin;tr.uLin(end)],'--','LineWidth',1);
                obj.decorate_time_axis(au,tc,refStart,distStart,'u [cm^3/s]','');

                % Osobny wykres pojedynczego scenariusza.
                f1 = obj.new_figure([100 100 960 620]);
                t1 = tiledlayout(f1,2,1,'TileSpacing','compact','Padding','compact');
                a1 = nexttile(t1); hold(a1,'on'); a1.Color = 'w';
                stairs(a1,tc,tr.yNL,'-','LineWidth',1);
                stairs(a1,tc,tr.yLin,'--','LineWidth',1);
                stairs(a1,tc,tr.ref,':','LineWidth',1);
                obj.decorate_time_axis(a1,tc,refStart,distStart,'y [cm]', ...
                    sprintf('Zadanie %+g cm, F_D %+g cm^3/s',scenarios(j,1),scenarios(j,2)));
                legend(a1,{'obiekt nieliniowy','obiekt liniowy','y_{zad}'},'Location','best');

                a2 = nexttile(t1); hold(a2,'on'); a2.Color = 'w';
                stairs(a2,tc,[tr.uNL;tr.uNL(end)],'-','LineWidth',1);
                stairs(a2,tc,[tr.uLin;tr.uLin(end)],'--','LineWidth',1);
                obj.decorate_time_axis(a2,tc,refStart,distStart,'u [cm^3/s]','');
                linkaxes([a1,a2],'x');
                obj.save_pdf(f1,sprintf('DMC_scenariusz_%d.pdf',j));
            end
            obj.save_pdf(fig,'regulacja_DMC_porownanie.pdf');
        end

        function plot_memberships(obj,x,weights,counts,tunedWeights)
            fig = obj.new_figure([100 100 980 680]);
            tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
            for m = 1:numel(counts)
                ax = nexttile(tiles); ax.Color = 'w';
                plot(ax,x,weights{m},'LineWidth',1);
                if counts(m)==5, title(ax,'5 reguł'); else, title(ax,sprintf('%d reguły',counts(m))); end
                xlabel(ax,'y [cm]'); ylabel(ax,'stopień przynależności \mu_i(y)');
                grid(ax,'on'); ylim(ax,[0 1.05]);
            end
            obj.save_pdf(fig,'funkcje_przynaleznosci.pdf');

            fig = obj.new_figure([100 100 760 420]);
            ax = axes(fig); ax.Color = 'w';
            plot(ax,x,tunedWeights,'LineWidth',1);
            title(ax,'Funkcje przynależności po strojeniu modelu TS (5 reguł)');
            xlabel(ax,'y [cm]'); ylabel(ax,'stopień przynależności \mu_i(y)');
            grid(ax,'on'); ylim(ax,[0 1.05]);
            obj.save_pdf(fig,'funkcje_przynaleznosci_strojone.pdf');
        end

        function plot_model_responses(obj,initial,tuned,steps,truth,t)
            selected = [1,6,7,12];
            fig = obj.new_figure([100 100 980 680]);
            tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
            for k = 1:numel(selected)
                j = selected(k);
                ax = nexttile(tiles); hold(ax,'on'); ax.Color = 'w';
                stairs(ax,t,truth(:,j),'LineWidth',1);
                stairs(ax,t,initial(:,j),'LineStyle','--','LineWidth',1);
                stairs(ax,t,tuned(:,j),'LineStyle','-.','LineWidth',1);
                grid(ax,'on'); xlabel(ax,'t [s]'); ylabel(ax,'y [cm]');
                title(ax,sprintf('Skok u %+g, F_D %+g cm^3/s',steps(j,1),steps(j,2)));
                if k==1
                    lg = legend(ax,{'obiekt nieliniowy','TS 5 początkowy','TS 5 po strojeniu'},'NumColumns',3);
                    lg.Layout.Tile = 'south';
                end
            end
            obj.save_pdf(fig,'porownanie_modeli_TS.pdf');
        end

        function plot_all_models(obj,linearOutputs,modelOutputs,counts,steps,truth,t)
            selected = [1,6,7,12];
            fig = obj.new_figure([100 100 1080 780]);
            tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
            for tile = 1:numel(selected)
                j = selected(tile);
                ax = nexttile(tiles); hold(ax,'on'); ax.Color = 'w';
                stairs(ax,t,truth(:,j),'LineWidth',1,'DisplayName','obiekt nieliniowy');
                stairs(ax,t,linearOutputs(:,j),'LineStyle',':','LineWidth',1,'DisplayName','liniowy');
                for m = 1:numel(counts)
                    stairs(ax,t,modelOutputs(:,j,m),'LineWidth',1, ...
                        'DisplayName',sprintf('TS %d',counts(m)));
                end
                grid(ax,'on'); xlabel(ax,'t [s]'); ylabel(ax,'y [cm]');
                title(ax,sprintf('Skok u %+g, F_D %+g cm^3/s',steps(j,1),steps(j,2)));
                if tile==1
                    lg = legend(ax,'NumColumns',6);
                    lg.Layout.Tile = 'south';
                end
            end
            obj.save_pdf(fig,'modele_2_5.pdf');
        end

        function plot_control_memberships(obj,x,initialWeights,tunedWeights,sigma)
            fig = obj.new_figure([100 100 960 390]);
            tiles = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
            titles = {'PDC: trójkątne funkcje przynależności jak w modelu', ...
                sprintf('Po strojeniu regulatora: funkcje dzwonowe (\\sigma = %.2f cm)',sigma)};
            data = {initialWeights, tunedWeights};
            for panel = 1:2
                ax = nexttile(tiles); ax.Color = 'w';
                plot(ax,x,data{panel},'LineWidth',1);
                title(ax,titles{panel});
                xlabel(ax,'y [cm]'); ylabel(ax,'stopień przynależności \mu_i(y)');
                grid(ax,'on'); ylim(ax,[0 1.05]);
            end
            obj.save_pdf(fig,'strojenie_regulatora.pdf');
        end

        function pdc(obj,t,traces,cases,refStart,distStart)
            labels = ["DMC cz. 1 (bez ogr.)"; "DMC ograniczony"; "PDC"; "PDC strojony"];
            styles = {'-','--','-','-.'};
            obj.plot_scenario_grid(t,traces,@(j) traces(j).y,labels,styles,cases, ...
                refStart,distStart,'y [cm]',true,'regulacja_porownanie.pdf',3);
            obj.plot_scenario_grid(t,traces,@(j) traces(j).u,labels,styles,cases, ...
                refStart,distStart,'u [cm^3/s]',false,'regulacja_sterowanie.pdf',4);
        end

        function plot_comparison(obj,t,traces,slTraces,cases,lambdaGrid,selectedVariant,refStart,distStart)
            slIdx = [2, selectedVariant, 4];
            labels = ["PDC"; "PDC strojony"; ...
                arrayfun(@(l) sprintf("FGPC-SL (\\lambda = %g)",l),lambdaGrid(slIdx).')];
            styles = {'--','-.',':','-','--'};
            obj.plot_scenario_grid(t,traces, ...
                @(j) [traces(j).y(:,3:4), slTraces(j,slIdx(1)).y, slTraces(j,slIdx(2)).y, slTraces(j,slIdx(3)).y], ...
                labels,styles,cases,refStart,distStart,'y [cm]',true,'porownanie_wyjsc.pdf',3);
            obj.plot_scenario_grid(t,traces, ...
                @(j) [traces(j).u(:,3:4), slTraces(j,slIdx(1)).u, slTraces(j,slIdx(2)).u, slTraces(j,slIdx(3)).u], ...
                labels,styles,cases,refStart,distStart,'u [cm^3/s]',false,'porownanie_sterowan.pdf',3);
        end

        function npl(obj,t,traces,slTraces,nplTraces,selectedVariant,lambdaGrid,cases,refStart,distStart)
            lamSel = lambdaGrid(selectedVariant);
            labels = ["DMC cz. 1"; "PDC strojony"; ...
                sprintf("FGPC-SL (\\lambda = %g)",lamSel); ...
                sprintf("FGPC-NPL (\\lambda = %g)",lamSel)];
            styles = {'--','-.','-','-'};
            obj.plot_scenario_grid(t,traces, ...
                @(j) [traces(j).y(:,1), traces(j).y(:,4), slTraces(j,selectedVariant).y, nplTraces(j,selectedVariant).y], ...
                labels,styles,cases,refStart,distStart,'y [cm]',true,'NPL_y.pdf',5);
            obj.plot_scenario_grid(t,traces, ...
                @(j) [traces(j).u(:,1), traces(j).u(:,4), slTraces(j,selectedVariant).u, nplTraces(j,selectedVariant).u], ...
                labels,styles,cases,refStart,distStart,'u [cm^3/s]',false,'NPL_u.pdf',4);
            obj.plot_tuning(t,traces,nplTraces,cases,lambdaGrid,refStart,distStart,'FGPC-NPL','strojenie_NPL.pdf');
        end

        function plot_tuning(obj,t,traces,fgpcTraces,cases,lambdaGrid,refStart,distStart,prefix,filename)
            % Wplyw parametru lambda dla duzego skoku dodatniego (+8 cm) i ujemnego (-8 cm).
            [~,jPos] = max(cases(:,1));
            [~,jNeg] = min(cases(:,1));
            picked = [jPos, jNeg];
            showIdx = [1, 3, 4];
            styles = {'--','-','-.'};
            labels = arrayfun(@(l) sprintf('%s (\\lambda = %g)',prefix,l),lambdaGrid(showIdx),'UniformOutput',false);

            fig = obj.new_figure([100 100 1060 720]);
            tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
            for col = 1:2
                j = picked(col);
                ay = nexttile(tiles,col); hold(ay,'on'); ay.Color = 'w';
                for idx = 1:numel(showIdx)
                    stairs(ay,t,fgpcTraces(j,showIdx(idx)).y,'LineStyle',styles{idx},'LineWidth',1);
                end
                stairs(ay,t,traces(j).ref,'LineStyle',':','LineWidth',1);
                obj.decorate_time_axis(ay,t,refStart,distStart,'y [cm]', ...
                    sprintf('Wyjście y (%s): zadanie %+g cm, F_D %+g cm^3/s',prefix,cases(j,1),cases(j,2)));

                au = nexttile(tiles,col+2); hold(au,'on'); au.Color = 'w';
                for idx = 1:numel(showIdx)
                    u = fgpcTraces(j,showIdx(idx)).u;
                    stairs(au,t,[u;u(end)],'LineStyle',styles{idx},'LineWidth',1);
                end
                obj.decorate_time_axis(au,t,refStart,distStart,'u [cm^3/s]', ...
                    sprintf('Sterowanie u (%s): zadanie %+g cm, F_D %+g cm^3/s',prefix,cases(j,1),cases(j,2)));

                if col==1
                    lg = legend(ay,[labels,{'y_{zad}'}],'NumColumns',4);
                    lg.Layout.Tile = 'south';
                end
            end
            obj.save_pdf(fig,filename);
        end
    end

    methods (Access = private)
        function plot_scenario_grid(obj,t,traces,curvesFn,labels,styles,cases,refStart,distStart,yLabel,withRef,filename,legendCols)
            fig = obj.new_figure([100 100 1100 840]);
            tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
            for j = 1:numel(traces)
                ax = nexttile(tiles); hold(ax,'on'); ax.Color = 'w';
                curves = curvesFn(j);
                if size(curves,1)==numel(t)-1
                    curves = [curves; curves(end,:)];
                end
                hL = stairs(ax,t,curves);
                for m = 1:numel(labels)
                    set(hL(m),'LineStyle',styles{m},'LineWidth',1);
                end
                if withRef
                    refL = stairs(ax,t,traces(j).ref,'LineStyle',':','LineWidth',1);
                end
                obj.decorate_time_axis(ax,t,refStart,distStart,yLabel, ...
                    sprintf('Zadanie %+g cm, F_D %+g cm^3/s',cases(j,1),cases(j,2)));
                if j==1
                    if withRef
                        lg = legend(ax,[hL(:);refL],[labels;"y_{zad}"],'NumColumns',legendCols);
                    else
                        lg = legend(ax,hL,labels,'NumColumns',legendCols);
                    end
                    lg.Layout.Tile = 'south';
                end
            end
            obj.save_pdf(fig,filename);
        end

        function save_pdf(obj,fig,filename)
            exportgraphics(fig,fullfile(obj.OutputDir,filename), ...
                'ContentType','vector','BackgroundColor','white');
            close(fig);
        end
    end

    methods (Static, Access = private)
        function fig = new_figure(pos)
            fig = figure('Visible','off','Color','w','Position',pos);
            theme(fig,'light');
        end

        function decorate_time_axis(ax,t,refStart,distStart,yLabel,titleText)
            xline(ax,refStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
            xline(ax,distStart,':','Color',[.45 .45 .45],'HandleVisibility','off');
            grid(ax,'on'); xlim(ax,[t(1),t(end)]);
            xlabel(ax,'t [s]'); ylabel(ax,yLabel);
            if ~isempty(titleText), title(ax,titleText); end
        end
    end
end
