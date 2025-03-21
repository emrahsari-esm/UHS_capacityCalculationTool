%% VE simulation in a standard black-oil solver
%  In this example we show how to set up a standard format black-oil
%  model that can be used to simulate a VE model. For the actual
%  simulation,  we use the fully-implicit solver in MRST from the 'ad-fi'
%  module, which is based on automatic differentiation.

try
   require deckformat ad-fi
catch %#ok<CTCH>
   mrstModule add deckformat ad-fi
end

%% Parameters for the simulation
gravity off
[nx,ny,nz] = deal(100, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(400,10,15); % Physical dimensions of reservoir
total_time = year/25;             % Total simulation time
nsteps     = 25;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 100;                % Permeability in milli darcies
phi        = 0.1;                % Porosity
depth      = 1000;               % Initial depth
ipress     = 200;                % Initial pressure

%% Create input deck and construct grid
% Create an input deck that can be used together with the fully-implicit
% solver from the 'ad-fi' module. Since the grid is constructed as part of
% setting up the input deck, we obtain it directly.
[deck, G] = sinusDeckAdi([nx ny nz], [Lx Ly H], nsteps, dt, ...
                         -.1*pi/180, depth, phi, perm, ...
                         (H*phi*Lx*Ly)*0.2*day/year, ipress);

% Alternatively, we could read deck from file and construct the grid
% deck = readEclipseDeck( ...
%    fullfile(VEROOTDIR,'data','decks','sinusDeckAdi.DATA');
% G = initEclipseGrid(deck);

figure, plotGrid(G),view([0 -1 0]), box on


%% Initialize data structures
% First, we convert the input deck to SI units, which is the unit system
% used by MRST. Second, we initialize the rock parameters from the deck;
% the resulting data structure may have to be post-processed to remove
% inactive cells. Then we set up the fluid object and tell the ad-fi solver
% that that we are working with an oil-gas system.
deck  = convertDeckUnits(deck);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);
% set the capillary pressure and the VE relperms explicitely
%Gt = topSurfaceGrid(G);
%fluid_case='hystersis';
fluid_case='simple';
switch fluid_case
    case 'simple'
       fluid.krW=@(sg,varargin) sg.^2;
       fluid.krOW=@(so,varargin) so.^2;
       fluid.krO = fluid.krOW;
       fluid.pcOW=@(sg,varargin) 1*barsa*(sg);
       mu1=0.4e-3;mu2=0.4e-4;
       fluid.muG=@(p,T) mu1+(mu2-mu1).*(T-273)./(300-273);
       mu1=1e-3;mu2=0.4e-3;
       fluid.muO=@(p,T) mu1+(mu2-mu1).*(T-273)./(300-273);
       fluid=rmfield(fluid,'relPerm');
       res_gas=0;
    otherwise
       disp('Use deck as fluid')
end


cW       = 4.1813*10.^6; % use volometric heatcapastity
cR       = 2.17*10.^6;
fluid.eW = @(T) cW.*T;
fluid.eO = @(T) cW.*T;
fluid.eR = @(T) cR.*T;
mycase='M1I2'
mycase='M1I1'
switch mycase
    case 'M1I1'
        % add reactions
        nm=1;
        ni=1;
        nr=1;
        % rate of right and left reaction
        fluid.LR=cell(nr,1);
        fluid.RR=cell(nr,1);
        % integers representing reaction equations
        fluid.ILn=zeros(nr,ni);
        fluid.IRn=zeros(nr,ni);
        fluid.Mn=zeros(nr,nm);
        % I1 <-> M1
        fluid.LR{1} = @(T)  exp(0.*T+log(1e-4));
        fluid.RR{1} = @(T)  exp((-273./T)-1+log(1e-4));
        fluid.ILn(1,1)=1;
        fluid.Mn(1,1)=1;
        M0=ones(G.cells.num,1).*0.0;
        M0(floor(G.cartDims(1)/2):floor(G.cartDims(1)*5/9),1)=0.5;
        I0=ones(G.cells.num,ni);
        I0(floor(G.cartDims(1)/2):floor(G.cartDims(1)*3/4),1)=0.1;
    case 'M1I2'
        nm=1;
        ni=2;
        % add reactions
        nr=1;
        % rate of right and left reaction
        fluid.LR=cell(nr,1);
        fluid.RR=cell(nr,1);
        % integers representing reaction equations
        fluid.ILn=zeros(nr,ni);
        fluid.IRn=zeros(nr,ni);
        fluid.IMn=zeros(nr,ni);
        % I1 <-> M1
        fluid.LR{1} = @(T)  exp(0.*T+log(1e-3));
        fluid.RR{1} = @(T)  exp((-273./T)-1+log(1e-11));
        fluid.ILn(1,1)=1;
        fluid.ILn(1,2)=2;
        fluid.Mn(1,1)=1;
        M0=ones(G.cells.num,1).*0.0;
        M0(floor(G.cartDims(1)/2):floor(G.cartDims(1)*5/9),1)=0.5;
        I0=zeros(G.cells.num,2);
        I0(floor(G.cartDims(1)/2):floor(G.cartDims(1)*3/4),1)=0.1;
        I0(floor(G.cartDims(1)/2):floor(G.cartDims(1)*4/5),2)=0.4;
      case 'M1I3'
        nm=1;
        ni=3;
        % add reactions
        nr=2;
        % rate of right and left reaction
        fluid.LR=cell(nr,1);
        fluid.RR=cell(nr,1);
        % integers representing reaction equations
        fluid.ILn=zeros(nr,ni);
        fluid.IRn=zeros(nr,ni);
        fluid.Mn=zeros(nr,nm);
        % ILn(1,1)*I1 + ILn(1,2)*I2 <-> M1
        fluid.LR{1} = @(T)  exp(0.*T+log(1e-3));
        fluid.RR{1} = @(T)  exp((-273./T)-1+log(1e-4));
        fluid.ILn(1,1)=1;
        fluid.ILn(1,2)=2;
        fluid.Mn(1,1)=1;
        % second reaction
        fluid.LR{2} = @(T)  exp(0.*T+log(1e-4));
        fluid.RR{2} = @(T)  exp((-273./T)-1+log(1e-4));
        fluid.ILn(2,1)=1;
        fluid.ILn(2,2)=1;
        fluid.IRn(2,3)=2;
        %fluid.Mn(1,1)=1;

        M0=ones(G.cells.num,1).*0.0;
        M0(floor(G.cartDims(1)/2):floor(G.cartDims(1)*5/9),1)=0.5;
        I0=zeros(G.cells.num,ni);
        I0(floor(G.cartDims(1)/2):floor(G.cartDims(1)*3/4),1)=0.1;
        I0(floor(G.cartDims(1)/2):floor(G.cartDims(1)*4/5),2)=0.4;
    otherwise
        error()
end


systemOW  = initADISystem({'Oil', 'Water','T','MI'}, G, rock, fluid);
% calculate conducivity for fock

fake_rock.perm=4.0*ones(G.cells.num,1);
T = computeTrans(G,fake_rock);
Trans=1./accumarray(G.cells.faces(:,1),1./T,[G.faces.num,1]);
internal=all(G.faces.neighbors>0,2);

systemOW.s.TH=Trans(internal);
systemOW.nonlinear.linesearch=true;

%% Run the schedule setup in the file
% Before we can run the schedule, we make sure that we have an initial
% hydrostatic pressure distribution. Then we pick the schedule from the
% input deck and start the simulator.
x0 = initEclipseState(G, deck, initEclipseFluid(deck));
z  = G.cells.centroids(:,3);
x0.pressure = ipress*barsa +(z(:)-z(end))*norm(gravity)*deck.PROPS.DENSITY(2);
x0.s(:,1)=deck.SOLUTION.SWAT;
x0.s(:,2)=deck.SOLUTION.SOIL;
x0.smax=x0.s;
x0.smin=x0.s;
x0.T=300*ones(G.cells.num,1);
x0.I=I0;
x0.M=M0;
%x0.T(floor(G.cartDims(1)/2))=330;

%x0.s=x0.s(:,[2,1]);
%x0.s=x0.z(:,[2,1]);

%Wext=processWells(G, rock, deck.SCHEDULE.control(1));
Wext=processWells(G, rock, deck.SCHEDULE.control(1));
for i=1:numel(Wext)
  Wext(i).T=273;
  Wext(i).I=0*ones(1,size(x0.I,2));
end
mrst_schedule = deck.SCHEDULE;
mrst_schedule.W={Wext};
%systemOW.fluid = fluid;
%systemOW.getEquations =@ eqsfiOWTMIExplicitWells;
%systemOW.stepFunction =@ stepOWTMI;
%[wellSols, states] = runScheduleADI(x0, G, rock, systemOW, deck.SCHEDULE,'Wext',Wext);
[wellSols, states] = runMrstADI(x0, G, systemOW, mrst_schedule,'Wext',Wext);
%[wellSols, states] = runScheduleADI(x0, G, rock, systemOW, deck.SCHEDULE,'Wext',Wext);

%% Plot results
%figure
figure(1)
xc = G.cells.centroids(:,1);
for nn=1:numel(states)
    clf
    state=states{nn};

        %
      figure(1)
      subplot(3,1,1),cla
      plot(xc,state.pressure/barsa);
      subplot(3,1,2),cla
      plot(xc,state.T);
      minT=min(state.T);maxT=max(state.T);
      subplot(3,1,3),cla
      plot(xc,state.s(:,2),xc,(state.T-minT)./(maxT-minT))
      hold on
      plot(xc,state.I,'x');
      plot(xc,state.M,'*');
      %figure(33)
      %plot(xc,state.M)
    drawnow;
    pause(1)
end

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
