classdef GasReservoirModel
    properties

        area;
        thickness;
        res_L;

        porosity;
        permeability;
        SWi;

        formation_top;
        rock_comp;

        res_temperature;

        initial_gas_composition;

        n_wells;
        wells_active;

        model_2d; % 2d model
        model_3d; % 3d model
        model_W;  % wells
        model_wellNo;
        model_rock;

        %output parameters
        model_pore_volume_total;
        model_z_factor_avg;
        model_ioip_ton;
        model_igip_sm3;
        model_p_avg;
        model_b0_avg;

        model_state;

    end

    methods
        % Constructor
        function obj = GasReservoirModel(area_, h_, porosity_, permeability_, SWi_, formation_top_, rock_comp_, n_wells_, res_temperature_, gas_composition_)

            % disp (area_);
            % disp (h_);
            % disp (porosity_);
            % disp (permeability_);
            % disp (SWi_);
            % disp (rock_comp_);
            % disp (n_wells_);
            % disp (res_temperature_);

            obj.area = area_;                   % in km2
            obj.thickness = h_;                 % in m
            obj.porosity = porosity_;       % in percent
            obj.permeability = permeability_;   % in mD
            obj.SWi = SWi_;                 % in percent
            obj.rock_comp = rock_comp_;         %in 1/bar
            obj.res_temperature = res_temperature_;     % in C

            %disp(gas_composition_);
            obj.initial_gas_composition = gas_composition_ / 100;
            %disp(obj.initial_gas_composition);

            obj.formation_top=formation_top_;   %in m

            obj.n_wells = n_wells_;

            obj.res_L = sqrt(obj.area*1e6);    %in m

            mrstModule add ad-core compositional deckformat ad-props mrst-gui

            % Dimensions
            nx = 1;
            ny = nx;
            nz = 10;

            % Set up grids
            dims = [nx, ny, nz];
            pdims = [obj.res_L, obj.res_L, obj.thickness];

            G = cartGrid(dims, pdims);

            % Reservoir geometry and petrophysical properties
            G    = computeGeometry(G);
            rock = makeRock(G, obj.permeability*milli*darcy, obj.porosity/100);

            obj.model_3d = G;
            obj.model_rock = rock;

        end

        % Function to build fluid model
        function obj=build_fluid_model(obj, gas_comp_, res_P_)

            gas_comp = gas_comp_;
            res_P = res_P_;

            % calculate_igip

            % initial settings
            sW_initial=obj.SWi/100;
            sO_initial = 0;
            sG_initial = 1-sW_initial-sO_initial;
            sW_irreducable=sW_initial;

            % Define fluid model
            gravity reset on

            % Surface Conditions
            p_sc = 101325*1; % atmospheric pressure
            T_sc = 293.15; % 20 C

            [mixture, info] = getBenchmarkMixture_custom( gas_comp, res_P*barsa, obj.res_temperature );
            %ncomp = mixture.getNumberOfComponents();

            eosname = 'pr'; % Select the Peng-Robinson equation of state
            EOSModel = EquationOfStateModel(obj.model_3d, mixture, eosname);

            %phase densities at surcafe cond.
            [~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

            nkr = 2;
            flowFluid = initSimpleADIFluid('phases', 'WOG', ...
                'n', [nkr, nkr, nkr], ...
                'rho', [1000, rhoO_S, rhoG_S], ...
                'mu', [1, 10, 0.05], ...
                'smin', [sW_irreducable,0.0,0.0],...
                'cR', obj.rock_comp/barsa,...
                'pRef' , p_sc);

            % 5: setup model

            % Construct models
            arg = {obj.model_3d, obj.model_rock, flowFluid, ...                   % Standard arguments
                mixture,...                                  % Compositional mixture
                'water', true, 'oil', true, 'gas', true};
            overallCompModel = GenericOverallCompositionModel(arg{:});   % Overall mole fractions

            % Validate model to initialize the necessary state function groups
            overallCompModel = overallCompModel.validateModel();

            % disp("model");

            if ~exist('maxChange', 'var')
                maxChange = 0.1;
            end
            if isfinite(maxChange)
                % Limit changes - natural variables limits phase mole fractions and
                % saturation while overall composition limits the overall mole fraction
                % change
                overallCompModel.dzMaxAbs = maxChange;
            end


            % Set up BC + initial conditions/initial guess
            % We start by the boundary conditions. Here, we use standard routines from
            % MRST, but must remember to also impose component specification at all
            % boundaries
            eos = overallCompModel.EOSModel;


            % We first set the initial state from the prescribed parameters, validate
            % the model (which will add data structures for any wells, etc), and then
            % perform a flash to ensure that the corresponding initial guess exists at
            % equilibrium conditions
            state0 = initCompositionalState(obj.model_3d, info.pressure, info.T, [sW_initial, sO_initial, sG_initial], info.initial, eos);
            state0 = overallCompModel.validateState(state0);

            % Set up bad initial guess
            initGuess = state0;
            %initGuess.components = zeros(1, eos.getNumberOfComponents());
            %initGuess.components(1) = 1;
            initGuess = overallCompModel.computeFlash(initGuess);

            % Adjust solver settings and construct nonlinear solver
            % To improve the nonlinear solver we impose limits on changes during the
            % nonlinear iterations. With natural variables, we limit phase mole
            % fractions and saturation, whereas for the overall composition we limit
            % the overall mole fraction
            % change
            if ~exist('maxChange', 'var')
                maxChange = 0.1;
            end
            if isfinite(maxChange)
                overallCompModel.dzMaxAbs = maxChange;
                %natural.dzMaxAbs = maxChange;
                %natural.dsMaxAbs = maxChange;
            end

            % Set up nonlinear solver
            nls = NonLinearSolver('reportLevel', 3,'MaxIterations', 1000);

            % Compute a single time step
            dt = 1000*day;
            %[state1_n, nreports]  = nls.solveTimestep(state0, dt, natural, ...
            %    'bc', [], 'initialGuess', initGuess);

            [state1_o, ~] = nls.solveTimestep(state0, dt, overallCompModel, ...
                'bc', [], 'initialGuess', initGuess);

            obj.model_state = state1_o;
            % IGIP

            [b0, pv] = overallCompModel.getProps(state1_o, 'ShrinkageFactors', 'PoreVolume');
            [rho, sat] = overallCompModel.getProps(state1_o, 'Density', 's');

            z_factor =  overallCompModel.getProps(state1_o, 'Z_V');

            sW0 = sat(:,1);
            sO0 = sat(:,2);
            sG0 = sat(:,3);
            pore_volume_total = sum(pv);

            obj.model_pore_volume_total =  pore_volume_total;

            obj.model_z_factor_avg = sum(z_factor.*pv)/pore_volume_total;

            obj.model_ioip_ton  = sum(pv.*(1-sW0-sG0).*rho{2})/1e3;

            obj.model_igip_sm3  = sum(pv.*sG0.*b0{3})/1e6;
            obj.model_p_avg = sum(state1_o.pressure.*pv)/pore_volume_total/barsa;

            obj.model_b0_avg = sum(b0{3}.*pv)/pore_volume_total;

            % disp("model_pore_volume_total");
            % disp(obj.model_pore_volume_total);

            % disp("model_z_factor_avg");
            % disp(obj.model_z_factor_avg);

            %disp("***************");
            %disp(state1_o);
            %disp("***************");

            %  disp(sum(state1_o.components(1)));
            %format long;
            %disp(state1_o.components(:,11));
            %disp(state1_o.pressure);

            % disp(pv);

            %disp("rock comp");
            %disp(obj.rock_comp);

        end

        function result = get_avg_res_pressure(obj)

            avg_p=obj.model_p_avg;
            result = avg_p;
        end

        function result = get_total_poreVolume(obj)

            pv=obj.model_pore_volume_total;
            result = pv;
        end

        function result = get_GIP(obj)

            gip=obj.model_igip_sm3;
            result = gip;
        end

        function result = get_z_factor(obj)

            z=obj.model_z_factor_avg;
            result = z;
        end

        function result = get_gas_Bg(obj)

            Bg=obj.model_b0_avg;
            result = Bg;
        end

        function result = get_initial_gas_composition(obj)

            x=obj.initial_gas_composition;
            result = x;
        end

        function plot_reservoir(obj, x_, y_, w_, h_, plot_no)

            width = w_;
            height = h_;
            window_x = x_;
            window_y = y_;

            margin_x =  width/4;
            margin_y = height/4;


            % Plotting functionality
            %fig2D = @() figure('Position', [(Pix_SS(3)-width)/2 + width/2 (Pix_SS(4)-height)/2 width height]);
            %fig3D = @() figure('Position', [(Pix_SS(3)-width)/2 + width (Pix_SS(4)-height)/2 width height]);

            fig2D = @() figure('Position', [window_x+margin_x window_y+margin_y width/2 height/2]);
            fig3D = @() figure('Position', [window_x+margin_x window_y+margin_y width/2 height/2]);


            alpha = 0.6;
            cmap  = jet*alpha + (1-alpha);
            setAxProps = @(ax) set(ax, 'View', [-45, 20], ...
                'PlotBoxAspectRatio', [2.00, 2.00, 1.00], ...
                'Projection'        , 'Perspective'     , ...
                'Box'               , 'on'              , ...
                'ColorMap'          , cmap              );

            pw = @() plotWell(obj.model_3d, obj.model_W, 'color', 'k', 'FontSize',9, 'height', 10, 'LineWidth', 2);

            if (plot_no == 1)  %2D

                fig2D(), plotGrid(obj.model_2d, 'EdgeAlpha', 0.1); axis equal tight, box on                 % Grid
                %plotFaces(obj.model_2d, obj.model_2d.faces.tag, 'edgeColor', 'b', 'lineWidth', 2); % Faults
                plotGrid(obj.model_2d, obj.model_2d.cells.tag , 'faceColor', 'r');

            elseif(plot_no==2) %3D
                fig3D(); plotGrid(obj.model_3d, 'EdgeAlpha', 0.1); %pw();
                setAxProps(gca), camlight();
                colormap(jet);
                title('3-D Plot of Reservoir [m]');
                axis on

            elseif(plot_no==3)  %pressure dist.

                fig3D(); plotCellData(obj.model_3d, obj.model_state.pressure/barsa, 'edgealpha', 0.2);
                %outlineCoarseGrid(G0, compartID,'EdgeColor','w','LineWidth',2);
                setAxProps(gca), camlight();
                colormap(jet);
                title('Pressure Distribution [bar]');
                axis on;
                colorbar;
            elseif(plot_no==4)  %h2 dist.

                fig3D(); plotCellData(obj.model_3d, obj.model_state.components(:, 11), 'edgealpha', 0.2);
                %outlineCoarseGrid(G0, compartID,'EdgeColor','w','LineWidth',2);
                setAxProps(gca), camlight();
                colormap(jet);
                title('Pressure Distribution [bar]');
                axis on;
                colorbar;
            end

        end

    end
end
