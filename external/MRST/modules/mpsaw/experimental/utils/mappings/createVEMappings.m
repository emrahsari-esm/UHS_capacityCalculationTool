function [cno, nno, hfno, fno,w]  = createVEMappings(g)
%Undocumented Utilty Function

%{ 
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
%}

   %% Create mapping from sub-half-face to cell, node, face, half-face and
   %% sub-face
   cellno   = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
   col      = 1+(cellno == g.faces.neighbors(g.cells.faces(:,1), 2));
   nhfaces  = g.cells.facePos(end)-1;
   hfaces   = accumarray([g.cells.faces(:,1), col], 1:nhfaces);
   hfaces   = rldecode(hfaces, diff(g.faces.nodePos));

   
   cells    =  rldecode(g.faces.neighbors, diff(g.faces.nodePos));
   nodes    =  repmat(g.faces.nodes, [2,1]);
   faces    =  repmat(rldecode(1:g.faces.num, diff(g.faces.nodePos),2)', [2,1]);
   %subfaces =  repmat((1:size(g.faces.nodes,1))', [2,1]);
   i        =  cells~=0;
   w        =  [cells(i), nodes(i), hfaces(i), faces(i)];% subfaces(i)];
   w        =  double(sortrows(w));
   
   
   cno     = w(:,1); 
   nno     = w(:,2);
   hfno    = w(:,3);
   fno     = w(:,4);
   %subfno  = w(:,5);
   %subhfno = (1:numel(cno))';
end
