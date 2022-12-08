function plot_Pareto(penalty_parameter,objective1,objective2,varargin)
% plot_Pareto plots the Pareto curves for a multi-objective optimization
% problem.
%
% ## Comments
%
% _none_
%
% ## Input Arguments
%
% `penalty_parameter` (_double_) - values of the penalty parameter
%
% `objective1` (_double_) - values of the first objective
%
% `objective2` (_double_) - values of the second objective
%
% `varargin` - optional arguments
%
% ## Output Arguments
%
% _none_
%

% default values
defaultNewfigure = true;
defaultLabels = {'Penalty Parameter' 'Objective 1' 'Objective 2'};

% input parser
p = inputParser;
addOptional(p,'newfigure',defaultNewfigure);
addOptional(p,'labels',defaultLabels);
parse(p,varargin{:});

% plot
if p.Results.newfigure
    figure
end
xlabel(p.Results.labels(1))
yyaxis left
plot(penalty_parameter,objective1);
ylabel(p.Results.labels(2))
yyaxis right
plot(penalty_parameter,objective2);
ylabel(p.Results.labels(3))

end