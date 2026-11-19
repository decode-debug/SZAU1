classdef ResultPlots
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
    end
end
