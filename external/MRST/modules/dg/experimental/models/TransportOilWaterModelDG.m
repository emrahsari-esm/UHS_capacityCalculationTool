classdef TransportOilWaterModelDG < TransportBlackOilModelDG
    % Two phase oil/water system without dissolution with discontinuous
    % Galerking discretization

    methods
        function model = TransportOilWaterModelDG(G, rock, fluid, varargin)
            model = model@TransportBlackOilModelDG(G, rock, fluid, varargin{:});
            model.oil   = true;
            model.water = true;
            model.gas   = false;
            model.conserveWater = true;
            model.conserveOil   = false;
        end

        % ----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] ...
                = transportEquationOilWaterDG(state0, state, model, dt, drivingForces, ...
                                  'solveForOil'  , model.conserveOil  , ...
                                  'solveForWater', model.conserveWater, ...
                                  varargin{:}                         );
        end
        
    end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}