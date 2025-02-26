classdef CoarseGrid1DCubature < Cubature
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = CoarseGrid1DCubature(G, prescision, internalConn)
            % Set up cubatureature
            
            % Basic properties handled by parent class
            cubature = cubature@Cubature(G, prescision, internalConn);
            cubature.dim       = 1;
            % Make cubature points and weights
            [x, w, pos] = cubature.makeCubature();
            % Assing properties
            cubature.points    = x;
            cubature.weights   = w;
            cubature.pos       = pos;
            
        end
        
        function [x, w, pos] = makeCubature(cubature)
           
            % Get cubature for parent grid
            parentCub = LineCubature(cubature.G.parent    , ...
                                     cubature.prescision  , ...
                                     cubature.internalConn);
            
            % Sort cubature points and weights according to partition
            
            order = cubature.G.faces.fconn;
            p     = rldecode((1:cubature.G.faces.num)', diff(cubature.G.faces.connPos), 1);
            [~, x, w]  = parentCub.getCubature(order, 'face');

            np  = diff(parentCub.pos);            
            np  = accumarray(p, np(order));
            pos = [0;cumsum(np)] + 1;
                        
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
