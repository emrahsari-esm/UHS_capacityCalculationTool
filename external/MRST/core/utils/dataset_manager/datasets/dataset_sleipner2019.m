function [info, present] = dataset_sleipner2019()
% Info function for SAIGUP dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

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

    helper = ['Visit the CO2DataShare website to download the dataset.',...
              ' To setup the full model grid, including feeder chimneys and wells,',...
              ' we require the dataset parts:',...
              ' Sleipner Reference Model 2019 Grid, Feeders, and Well data.'...  
              ' Extract the contents of the folders from the zip archive'...
              ' into a folder named "Sleipner2019" under the current data directory: "',...
              mrstDataDirectory(), '"'];

[info, present] = datasetInfoStruct(...
        'name', 'Sleipner2019', ...
        'website', 'https://co2datashare.org/dataset/sleipner-2019-benchmark-model', ...
        'fileurl', '', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', false, ...
        'cells', 1986176, ...
        'instructions', helper , ...        
        'examples', { ...
                    }, ...
    'description', ['The Sleipner 2019 Benchmark Model is a reference dataset from the Sleipner CO2 storage site. This site, located off the western coast of Norway using part of the larger Utsira Formation as a storage unit, has been actively used by Equinor as a CO2 storage site since 1996, for the purpose of avoiding emissions associated with the production of natural gas in the same region. The site has been extensively monitored from the beginning, and a series of time-lapse seismic datasets that documents the migration of the injected CO2 has been established and (along with associated well logs and baseline seismic) has served as input to the creation of the benchmark model.',...
    'The dataset includes a simulation grid with associated petrophysical properties and well data, as well as velocity maps and individual horizons within the storage reservoir. It was prepared over the course of 2018 and 2019, and officially released January 2020 under the SLEIPNER CO2 REFERENCE DATASET LICENSE, for the purpose of supporting research and development into CO2 storage issues. It succeeds an earlier reference model, the 2011 Sleipner Benchmark Model which over the years has been used by over 40 research institutes and served as a basis for a large amount of published research.'], ...
        'modelType', 'Corner-point' ...
         );
end
