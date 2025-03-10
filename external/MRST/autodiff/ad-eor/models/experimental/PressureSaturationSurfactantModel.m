classdef PressureSaturationSurfactantModel < OilWaterSurfactantBaseModel
%
%
% SYNOPSIS:
%   model = PressureSaturationSurfactantModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Model used to solve the pressure-saturation equation implicitly
% for a constant (in time) surfactant concentration.
%
% PARAMETERS:
%   G        - Grid
%   rock     - Rock structure
%   fluid    - Fluid structure
%   varargin - optional parameter
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO: ImplicitExplicitOilWaterSurfactantModel, equationsPressureSaturationForOilWaterSurfactant
%


    methods
        function model = PressureSaturationSurfactantModel(G, rock, fluid, varargin)
            model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsPressureSaturationForOilWaterSurfactant(state0, state, model, ...
                                                              dt, drivingForces, varargin{:});
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
