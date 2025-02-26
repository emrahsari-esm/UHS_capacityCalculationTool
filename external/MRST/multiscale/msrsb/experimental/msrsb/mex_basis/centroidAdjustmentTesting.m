mrstModule add multiscale-devel coarsegrid 

if 0
    G = cartGrid([60 5]);
    p = partitionUI(G, [3 1]);
    refined = p == 1;
    p(refined) = (1:sum(refined))' + 9;
    targetBlock = 1;
else
    G = cartGrid([60 60]);
    p = partitionUI(G, [3 3]);
    refined = p == 4;
    p(refined) = (1:sum(refined))' + 9;
    
    targetBlock = 4;
end


close all
plotToolbar(G, mod(p, 13))

G = computeGeometry(G);

rock.perm = repmat(500*milli*darcy, G.cells.num, 1);
A = getIncomp1PhMatrix(G, computeTrans(G, rock));




p = processPartition(G, p);
p = compressPartition(p);

CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG);

%%
sg = setupGridsForMex(A, CG);
active = true(G.cells.num, 1);

% active = refined;


msdir = mrstPath('query', 'multiscale-devel');
delete(fullfile(msdir, 'mex_basis', 'mex_iteratedJacobiBasisFaster.mex*'))
rehash

[I, I_comp] = mbasisSmoothed(sg, 'tolerance',2, 'maxiter', 1000, 'omega', .66, 'active', active);


close all
plotToolbar(G, I(:, targetBlock));
axis equal tight off
outlineCoarseGrid(G, p)
plotGrid(G, CG.cells.centers)
% clear tmp
% tmp.sum = sum(I, 2);
% tmp.I = I;

% close all
% plotToolbar(G, tmp)
% colorbar
%%
% sg = setupGridsForMex(A, CG);

close all
for i = 1:CG.cells.num
    figure(1); clf;
    
    c = sg.cells{i} + 1;
    cbnd = c(sg.isBnd{i} == 1);
    
    c = CG.cells.interaction{i};
    
    plotGrid(G, c, 'FaceColor', 'red');
%     plotGrid(G, setdiff(c, cbnd), 'FaceColor', 'red');
    plotGrid(G, cbnd, 'FaceColor', 'blue', 'FaceAlpha', .2);
    outlineCoarseGrid(G, p);
    pause()
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
