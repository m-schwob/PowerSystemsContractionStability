classdef DisplayTool
    methods (Static)
        % plot formula formatted as sym but add pre to latex
        % use: display_symbolic_pre({f1,'pre1'},{f2,'pre2'},...)
        function display_symbolic_pre(varargin)
            formulas = cell([1,nargin]);
            for k = 1:nargin
                formulas{k} = strcat(varargin{k}{2}, latex(varargin{k}{1}));       
            end
            DisplayTool.display_latex(formulas{:});
        end

        % plot formula formatted as sym
        function display_symbolic(varargin)
            formulas = cell([1,nargin]);
            for k = 1:nargin
                formulas{k} = latex(varargin{k});       
            end
            DisplayTool.display_latex(formulas{:});
        end

        % plot formulas formatted as latex
        function display_latex(varargin)
            figure; hold on; hold off;
            for k = 1:nargin
                subplot(nargin,1,k); hold all;
                axis off;
                text(0.5, 0.5, ['$$' varargin{k} '$$'], 'Interpreter','latex', 'FontSize',15, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle');
            end
            grid off; hold off;
        end
    end
end