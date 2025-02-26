function problemCopy = copyPackedProblem(problem, newBaseName, varargin)
% Make a copy of a packed problem with a new name and updated
% resultHandlers.
%
% SYNOPSIS:
%   problemCopy = copyPackedProblem(problem, newBaseName, varargin)
%
% REQUIRED PARAMETERS:
%   problem      - Existing problem to be copied. Generated by
%                  simulatePackedProblem.
%
%   newBaseName  - Updated BaseName for the new problem. This should be 
%                  different from the BaseName of the original problem
%                  unless a different 'Name' is given as an optional
%                  argument.
%
%
% OPTIONAL PARAMETERS:
%
%   Directory       - Folder where simulation output is to be stored. The
%                     default is fullfile(mrstOutputDirectory(), BaseName).
%
%   Name            - Name of this specific problem. We assume that
%                     BaseName + Name uniquely identifies all inputs to
%                     this function. Output will be stored under
%                     BaseName/Name.
%
%   Description     - A description of the case. Can be useful to add
%                     additional information to the shorted Name +
%                     BaseName pair.
%
%
% RETURNS:
%   problemCopy - A copy of the original problem with updated identifier
%                 (BaseName + Name) and updated ResultHandler
%                 dataDirectory.
%
% SEE ALSO:
%   packSimulationProblem, getPackedSimulatorOutput, simulatePackedProblem

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

    opt = struct('Directory',          '', ...
                    'Name',            '',...
                    'Description',     '');
             
    opt = merge_options(opt, varargin{:});
    
    % Check for duplicate names
    
    noName = isempty(opt.Name);
    if noName
        opt.Name = problem.Name;
    end

    if isempty(opt.Description)
        opt.Description = problem.Description;
    end

    if isempty(opt.Directory)
        opt.Directory = fullfile(mrstOutputDirectory(), newBaseName);
    end
    
    
    if strcmp(problem.BaseName,newBaseName) && strcmp(problem.Name,opt.Name)
        warning('Original problem and copy have the same name.')
        if strcmp(problem.OutputHandlers.states.getDataPath,fullfile(opt.Directory,opt.Name))
            error('Cannot create new problem with same output directory as original problem.')
        end
    end
    
    problemCopy = problem;
    
    problemCopy.BaseName = newBaseName;
    problemCopy.Name = opt.Name;
    problemCopy.Description = opt.Description;
    
    
    makeHandler = @(prefix) ResultHandler('dataPrefix', prefix, ...
                                          'writeToDisk', true,...
                                          'dataDirectory', opt.Directory, ...
                                          'dataFolder', problemCopy.Name, ...
                                          'cleardir', false);

    problemCopy.OutputHandlers.states   = makeHandler('state');
    problemCopy.OutputHandlers.reports  = makeHandler('report');
    problemCopy.OutputHandlers.wellSols = makeHandler('wellSols');
    
end