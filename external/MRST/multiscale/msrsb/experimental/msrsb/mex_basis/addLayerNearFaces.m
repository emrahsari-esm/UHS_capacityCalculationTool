function [p, cells] = addLayerNearFaces(G, p, faces, additive)
    % Add extra partitions near certain fine faces

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


    if nargin < 4
        additive = true;
    end
    isBF = false(G.faces.num, 1);
    isBF(boundaryFaces(G)) = true;
    
    assert(all(isBF(faces)));
    
    cells = sum(G.faces.neighbors(faces, :), 2);
    
    if additive
        p(cells) = max(p) + p(cells);
    else
        p(cells) = max(p) + 1;
    end
    p = compressPartition(p);
end
