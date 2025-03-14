function mrsttovtk(G, states, name, formatn)
%  Write one or more data sets into the files name.pvd and name-n.vtu 
%  for visualization in ParaView.
%
% SYNOPSIS:
%       mrsttovtk(G, struct('data1', data, 'data2', data), name, formatn)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   data    - Scalar cell data.
%
%   name    - Name of the outputs files name.pvd and name-n.vtu.
%
%   formatn - Format for writing numerical values.
%
% NOTES:
%   The function `mrsttovtk` only writes cell values (e.g., pressures) but 
%   currenlty not values on faces (e.g., fluxes). This function supports 3D 
%   tetrahedral, polyhedral, and corner-pont grids, 2D cartesian, triangular,
%   and polygonal grids, and 1D grids.
%
% EXAMPLE:
%   Let us consider the grid from the MRST 2019 book (pp 97, Fig 3.31):
%   G = simpleGrdecl([20, 20, 5], @(x) .25*(x-.5),'flat',true);
%   G = processGRDECL(G);
%
%   % 1) Writting the content in G.cells.indexMap into the files states.pvd
%   and states-00001.vtu for visualization in ParaView (this since indexMap
%   is the only field with the same numbers of values as grid cells).
%   mrsttovtk(G,G.cells,'states','%.6g');
%
%   % 2) After making a structure with 100 values, we run the mrsttovtk 
%   function to write the values to be visualized each 10 values in ParaView.
%   for i=1:100
%       states{i}.index=i^2*G.cells.indexMap;
%   end
%   mrsttovtk(G,{states{10:10:100}},'states','%.6g');

%{
Copyright 2021-2022, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling. 
This file is part of the ad-micp module.
ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
    clear Text
    form = sprintf('%s',formatn);
    forms = sprintf('%s ',formatn);
    Text{1} = sprintf('\t\t\t<Points>\n');
    Text{size(Text,2)+1} = sprintf(['\t\t\t\t<DataArray type="Float32" ' ...
            'Name="Coordinates" NumberOfComponents="3" format="ascii">\n']);
    Text{size(Text,2)+1}=sprintf('\t\t\t\t\t');
    if G.griddim == 1
        tg = 3;
        for i = 1 : G.nodes.num - 1
            Text{size(Text, 2) + 1} = sprintf([forms forms forms], ...
              G.nodes.coords(i, 1), 0, 0);
        end
        Text{size(Text,2)+1}=sprintf([forms forms form '\n'], ...
        G.nodes.coords(end,1),0,0);
    elseif G.griddim == 2
        if 4*G.cells.num==size(G.cells.faces, 1)
            tg = 8;
        else
            tg = 7;
        end
        for i=1:G.nodes.num-1
            Text{size(Text,2)+1}=sprintf([forms forms forms], ...
              G.nodes.coords(i,1),G.nodes.coords(i,2),0);
        end
        Text{size(Text,2)+1}=sprintf([forms forms form '\n'], ...
        G.nodes.coords(end,1),G.nodes.coords(end,2),0);
    else
        tg = 42;
        for i=1:G.nodes.num-1
            Text{size(Text,2)+1}=sprintf([forms forms forms], ...
              G.nodes.coords(i,1),G.nodes.coords(i,2),G.nodes.coords(i,3));
        end
        Text{size(Text,2)+1}=sprintf([forms forms form '\n'], ...
        G.nodes.coords(end,1),G.nodes.coords(end,2),G.nodes.coords(end,3));
    end
    Text{size(Text,2)+1}=sprintf('\t\t\t\t</DataArray>\n');
    Text{size(Text,2)+1}=sprintf('\t\t\t</Points>\n');
    Text{size(Text,2)+1}=sprintf('\t\t\t<Cells>\n');
    Text{size(Text,2)+1}=sprintf(['\t\t\t\t<DataArray type="UInt8" ' ...
                 'Name="types" NumberOfComponents="1" format="ascii">\n']);
    Text{size(Text,2)+1}=sprintf('\t\t\t\t\t');
    Text{size(Text,2)+1}=sprintf(repmat([num2str(tg),' '], 1, G.cells.num-1));
    Text{size(Text,2)+1}=sprintf('%d\n',tg);
    Text{size(Text,2)+1}=sprintf('\t\t\t\t</DataArray>\n');
    ii=rldecode(1:G.cells.num, diff(G.cells.facePos), 2)==1;
    if G.griddim == 3
        Text{size(Text,2)+1}=sprintf(['\t\t\t\t<DataArray type="Int32" ' ...
                 'Name="faces" NumberOfComponents="1" format="ascii">\n']);
        Text{size(Text,2)+1}=sprintf('\t\t\t\t\t');
        Text{size(Text,2)+1}=sprintf('%d',sum(ii));
        kk=-1;
        mm=0;
        for j=G.cells.faces(ii,1)'
            jj=rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)==j;
            Text{size(Text,2)+1}=sprintf(repmat(' %d', 1,...
                           sum(jj)+1),sum(jj),G.faces.nodes(jj,1)-1);
            ll=unique([G.faces.nodes(jj,1)'-1 kk]);
            kk=ll;
            mm=sum(jj)+mm;
        end 
        offsets{1}=ll(2:end);
        faceoffsets{1}=mm+sum(ii)+1;
        for i=2:G.cells.num
            ii=rldecode(1:G.cells.num, diff(G.cells.facePos), 2)==i;
            Text{size(Text,2)+1}=sprintf(' %d',sum(ii));
            kk=-1;
            mm=0;
            for j=G.cells.faces(ii,1)'
                jj=rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)==j;
                Text{size(Text,2)+1}=sprintf(repmat(' %d', 1,...
                           sum(jj)+1),sum(jj),G.faces.nodes(jj,1)-1);
                ll=unique([G.faces.nodes(jj,1)'-1 kk]);
                kk=ll;
                mm=sum(jj)+mm;
            end 
            offsets{i}=ll(2:end);
            faceoffsets{i}=mm+sum(ii)+1;
        end
        Text{size(Text,2)+1}=sprintf('\n');
        Text{size(Text,2)+1}=sprintf('\t\t\t\t</DataArray>\n');
    else
        kk=-1;
        mm=0;
        vct=zeros(max(size(G.cells.faces(ii,1))),1);
        jjj=G.cells.faces(ii,1)';
        jj=rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)==jjj(1);
        [ll,id,ic]=unique([G.faces.nodes(jj,1)'-1 kk],"stable");
        c = setdiff (ll, kk);
        vct(1)=c(1);
        vct(2)=c(2);
        kk=ll;
        mm=sum(jj)+mm;
        itr=3;
        for indx=2:size(jjj, 2)-1
            j=jjj(indx);
            jj=rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)==j;
            [ll,id,ic]=unique([G.faces.nodes(jj,1)'-1 kk],"stable");
            c = setdiff (ll, kk);
            vct(itr)=c;
            kk=ll;
            mm=sum(jj)+mm;
            itr=itr+1;
        end 
        offsets{1}=vct;
        faceoffsets{1}=mm+sum(ii)+1;
        for i=2:G.cells.num
            ii=rldecode(1:G.cells.num, diff(G.cells.facePos), 2)==i;
            kk=-1;
            mm=0;
            vct=zeros(max(size(G.cells.faces(ii,1))),1);
            jjj=G.cells.faces(ii,1)';
            jj=rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)==jjj(1);
            [ll,id,ic]=unique([G.faces.nodes(jj,1)'-1 kk],"stable");
            c = setdiff (ll, kk);
            vct(1)=c(1);
            vct(2)=c(2);
            kk=ll;
            mm=sum(jj)+mm;
            itr=3;
            for indx=2:size(jjj, 2)-1
                j=jjj(indx);
                jj=rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)==j;
                ll=unique([G.faces.nodes(jj,1)'-1 kk],"stable");
                c = setdiff (ll, kk);
                vct(itr)=c;
                kk=ll;
                mm=sum(jj)+mm;
                itr=itr+1;
            end 
            offsets{i}=vct;
            faceoffsets{i}=mm+sum(ii)+1;
        end
    end
    Text{size(Text,2)+1}=sprintf(['\t\t\t\t<DataArray type="Int32" ' ...
               'Name="offsets" NumberOfComponents="1" format="ascii">\n']);
    ii=0;
    Text{size(Text,2)+1}=sprintf('\t\t\t\t\t');
    for i=1:G.cells.num-1
        ii=max(size(offsets{i}))+ii;
        Text{size(Text,2)+1}=sprintf('%d ',ii);    
    end 
    Text{size(Text,2)+1}=sprintf('%d\n',max(size(offsets{G.cells.num})) ...
                                                                      +ii);  
    Text{size(Text,2)+1}=sprintf('\t\t\t\t</DataArray>\n');
    Text{size(Text,2)+1}=sprintf(['\t\t\t\t<DataArray type="Int32" '...
          'Name="connectivity" NumberOfComponents="1" format="ascii">\n']);
    Text{size(Text,2)+1}=sprintf('\t\t\t\t\t');
    for i=1:G.cells.num-1
        Text{size(Text,2)+1}=sprintf(repmat('%d ', 1, ...
                                  max(size(offsets{i}))),offsets{i});   
    end
    Text{size(Text,2)+1}=sprintf(repmat('%d ', 1, ...
        max(size(offsets{G.cells.num}))-1),offsets{G.cells.num}(1:end-1));
    Text{size(Text,2)+1}=sprintf('%d\n',offsets{G.cells.num}(end));                   
    Text{size(Text,2)+1}=sprintf('\t\t\t\t</DataArray>\n');
    if G.griddim == 3
        Text{size(Text,2)+1}=sprintf(['\t\t\t\t<DataArray type="Int32" '...
           'Name="faceoffsets" NumberOfComponents="1" format="ascii">\n']);
	    Text{size(Text,2)+1}=sprintf('\t\t\t\t\t');
        ii=0;
        for i=1:G.cells.num-1
            ii=faceoffsets{i}+ii;
            Text{size(Text,2)+1}=sprintf('%d ',ii);
        end
        Text{size(Text,2)+1}=sprintf('%d\n',faceoffsets{G.cells.num}+ii);
        Text{size(Text,2)+1}=sprintf('\t\t\t\t</DataArray>\n');
    end
    Text{size(Text,2)+1}=sprintf('\t\t\t</Cells>\n');
    Text{size(Text,2)+1}=sprintf('\t\t</Piece>\n');
    Text{size(Text,2)+1}=sprintf('\t</UnstructuredGrid>\n');
    Text{size(Text,2)+1}=sprintf('</VTKFile>\n');
    Text = [Text{:}];
    if max(size(states))>1
        state0=states{1};
    else
        state0=states;
    end
    fields=fieldnames(state0);
    CellData=fields;
    for m=1:max(size(fields))
        if max(size(getfield(state0,fields{m})))~=G.cells.num
            CellData(ismember(CellData,fields{m}))=[];
        end
    end
    filename = sprintf('%s.pvd',name);
    fileIDpvd = fopen(filename,'w');
    fprintf(fileIDpvd,'<?xml version="1.0"?>\n');
    fprintf(fileIDpvd,'<VTKFile type="Collection"\n');
    fprintf(fileIDpvd,'\t\t\t\t version="0.1"\n');
    fprintf(fileIDpvd,'\t\t\t\t byte_order="LittleEndian"\n');
    fprintf(fileIDpvd,'\t\t\t\t compressor="vtkZLibDataCompressor">\n');
    fprintf(fileIDpvd,'\t<Collection>\n');
    
    filename = sprintf('%s-%.05d.vtu',name,1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'<?xml version="1.0"?>\n');
    fprintf(fileID,['<VTKFile type="UnstructuredGrid" version="0.1" '...
                                          'byte_order="LittleEndian">\n']);
    fprintf(fileID,'\t<UnstructuredGrid>\n');
    fprintf(fileID,['\t\t<Piece NumberOfCells="%d" '...
                        'NumberOfPoints="%d">\n'],G.cells.num,G.nodes.num);
    fprintf(fileID,'\t\t\t<CellData Scalars="Variables">\n');
    for m=1:size(CellData,1)
        values=getfield(state0,CellData{m});
        if min(size(values))~=1
            for o=1:min(size(values))
                fprintf(fileID,['\t\t\t\t<DataArray type="Float32" '...
                'Name="%s-%d" NumberOfComponents="1" format="ascii">\n']...
                                                           ,CellData{m},o);
                fprintf(fileID,'\t\t\t\t\t');                                       
                fprintf(fileID,repmat(forms,1,G.cells.num-1), ...
                                                        values(1:end-1,o));
                fprintf(fileID,[form '\n'],values(end,o));
                fprintf(fileID,'\t\t\t\t</DataArray>\n');
            end
        else
            fprintf(fileID,['\t\t\t\t<DataArray type="Float32" '...
            'Name="%s" NumberOfComponents="1" format="ascii">\n'], ...
                                                              CellData{m});
            fprintf(fileID,'\t\t\t\t\t');
            fprintf(fileID,repmat(forms, 1,G.cells.num-1), ...
                                                        values(1:end-1,1));
            fprintf(fileID,[form '\n'],values(end,1));
            fprintf(fileID,'\t\t\t\t</DataArray>\n');
        end
    end
    fprintf(fileID,'\t\t\t</CellData>\n');
    fprintf(fileID,Text);
    fclose(fileID);
    if isfield(state0,'time')==1
        fprintf(fileIDpvd,'\t\t<DataSet timestep="%d" file="%s"/>\n', ...
                                                     state0.time,filename);
    else
        fprintf(fileIDpvd,'\t\t<DataSet timestep="0" file="%s"/>\n', ...
                                                                 filename);
    end
    for n=2:max(size(states))
        filename = sprintf('%s-%.05d.vtu',name,n);
        fileID = fopen(filename,'w');
        fprintf(fileID,'<?xml version="1.0"?>\n');
        fprintf(fileID,['<VTKFile type="UnstructuredGrid" '...
                            'version="0.1" byte_order="LittleEndian">\n']);
        fprintf(fileID,'\t<UnstructuredGrid>\n');
        fprintf(fileID,['\t\t<Piece NumberOfCells="%d" '...
                        'NumberOfPoints="%d">\n'],G.cells.num,G.nodes.num);
        fprintf(fileID,'\t\t\t<CellData Scalars="Variables">\n');
        for m=1:size(CellData,1)
            values=getfield(states{n},CellData{m});
            if min(size(values))~=1
                for o=1:min(size(values))
                    fprintf(fileID,['\t\t\t\t<DataArray type="Float32" '...
                    'Name="%s-%d" NumberOfComponents="1" '...
                                       'format="ascii">\n'],CellData{m},o);
                    fprintf(fileID,'\t\t\t\t\t');
                    fprintf(fileID,repmat(forms, 1, G.cells.num-1), ...
                                                        values(1:end-1,o));
                    fprintf(fileID,[form '\n'],values(end,o));
                    fprintf(fileID,'\t\t\t\t</DataArray>\n');
                end
            else
                fprintf(fileID,['\t\t\t\t<DataArray type="Float32" ' ...
                 'Name="%s" NumberOfComponents="1" format="ascii">\n'], ...
                                                              CellData{m});
                fprintf(fileID,'\t\t\t\t\t');                                          
                fprintf(fileID,repmat(forms, 1,G.cells.num-1), ...
                                                        values(1:end-1,1));
                fprintf(fileID,[form '\n'],values(end,1));
                fprintf(fileID,'\t\t\t\t</DataArray>\n');
            end
        end
        fprintf(fileID,'\t\t\t</CellData>\n');
        fprintf(fileID,Text);
        fclose(fileID);
        if isfield(state0,'time')==1
            fprintf(fileIDpvd,['\t\t<DataSet timestep="%d" '...
                                 'file="%s"/>\n'],states{n}.time,filename);
        else
            fprintf(fileIDpvd,['\t\t<DataSet timestep="%d" '...
                                            'file="%s"/>\n'],n-1,filename);
        end
    end
    fprintf(fileIDpvd,'\t</Collection>\n');
    fprintf(fileIDpvd,'</VTKFile>\n');
    fclose(fileIDpvd);
end
