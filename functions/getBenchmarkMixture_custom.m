function [fluid, info] = getBenchmarkMixture_custom( initial_composition, res_initial_pressure, res_initial_temp, varargin)

        names = {'Methane', 'Ethane', 'n-Propane', 'IsoButane', 'n-Butane', 'IsoPentane', 'n-Pentane', 'n-Hexane', 'CarbonDioxide', 'Nitrogen', 'Hydrogen' };

        info = makeInfo('injection', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0], ...
            'initial',   initial_composition, ...
            'pressure', res_initial_pressure, ...
            'T', 273.15 + res_initial_temp);

        fluid = TableCompositionalMixture(names);

%fluid.name = "descr";
end

function info = makeInfo(varargin)
info = struct('injection', [], 'initial', [],...
    'pressure',  [], 'temp',    [], 'T', []);
info = merge_options(info, varargin{:});
info.temp = info.T;
end