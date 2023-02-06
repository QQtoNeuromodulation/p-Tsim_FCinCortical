close all;
clear;
clc;
delete input.h5;
delete output.h5;
delete start_index.txt;
delete pre.mat;
delete Transducer_PosN.vtk;
delete Pressure_PosN.vtk;
delete PosN_acoustic.mat;
delete preThermal.mat;
delete T_max_PosN.vtk;
delete Q_PosN.vtk;
delete PosN_thermal.mat;

resolution = 0.3;
dt = 27e-9;
t_end = 75e-6;    % total compute time [s] (this must be long enough to reach steady state)

% =========================================================================
% Load CT image
% =========================================================================
Image_Original = niftiread('../MRI_CT.nii.gz');
%Image_Original = squeeze(Image_Original);
%Image_higres = Image_Original(211:512,111:408,91:325);
%Image_resized = imresize3(Image_higres,roundEven([0.234/resolution 0.234/resolution 0.2/resolution].*(size(Image_higres))));
%Image = permute(Image_resized,[1 3 2]);
Image = permute(Image_Original,[2 3 1]);
Image = flip(Image,2);
Image = flip(Image,3);
%X_Start = 100;
%X_End = 317;
%X_Expand = 110;
%Z_Start = 1;
%Z_Expand = 60;
%Y_Expand = 95;
%Y_End = 285;
%Image = expandMatrix(Image, [0, X_Expand, Y_Expand, 0, 0, Z_Expand]);
%Image = Image(X_Start:end,1:Y_End,Z_Start:end);

% Focus and Bowl coordinates determined in Slicer 
Focus_real.x = 249; 
Focus_real.y = 204;
Focus_real.z = 59;
Bowl_real.x = 150;
Bowl_real.y = 204;
Bowl_real.z = 248;

Focus.x = Focus_real.y + 1; % +1 added for MATLAB
Focus.y = length(Image_Original(1,1,:))-(Focus_real.z+1)+1;
Focus.z = length(Image_Original(:,1,1))-(Focus_real.x+1)+1;

Bowl.x = Bowl_real.y + 1; % +1 added for MATLAB
Bowl.y = length(Image_Original(1,1,:))-(Bowl_real.z+1)+1;
Bowl.z = length(Image_Original(:,1,1))-(Bowl_real.x+1)+1;
% To visualize in Paraview, use the converted Focus and Bowl and subtract 1
Picked_Slice = Focus.x;

% %% going backward
% X = 161.92541063011024;
% Y = 83.60365510218355;
% Z = 227.05175109498379;
% XX = length(Image_Original(:,1,1)) - Z -(Z_Start-1)
% YY = X + X_Start - 2
% ZZ = length(Image_Original(1,1,:)) - Y + Y_Expand

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% source parameters
source_f0       = 500e3;      % source frequency [Hz]
source_roc      = 63.2e-3;    % bowl radius of curvature [m]
source_diameter = 64e-3;    % bowl aperture diameter [m]
source_amp      = 37.288e3;      % source pressure [Pa]


% medium parameters
c_water         = 1482.3;     % sound speed [m/s]
c_skin          = 1624.0;
c_muscle        = 1588.4;
c_skull_max     = 3360.0;   % Skull cortical (max)
c_skull_min     = 1854.7;   % Skull Cancellous (min)

rho_water       = 998.23;     % density [kg/m^3]
rho_skin        = 1109;
rho_muscle      = 1090;
rho_skull_max   = 2100;
rho_skull_min   = 1080;

alpha_power_skull     = 1.1;          % b: water: 2      tissue: 1.1

% (alpha[Np/m]*8.686) / ((1/(1/Frequency[MHz]))*100) Converting Np/m to
% [dB/MHz/cm]

alpha_coeff_water     = 0.0022;
alpha_coeff_skull     = 8.0;
alpha_coeff_skin      = 1.84;
alpha_coeff_muscle    = 0.59;

% computational parameters
record_periods  = 1;        % number of periods to record

% =========================================================================
% RUN SIMULATION
% =========================================================================

% calculate the grid spacing based on the PPW and F0
dx = resolution*1e-3;

% compute the size of the grid
[Nx, Ny, Nz] = size(Image);

% --------------------
% Skull
% --------------------
% segment the skull based on the Hounsfield units
skull_mask = logical(flip(flip(permute(niftiread('../CT_Skull_mask.nii.gz'),[2 3 1]),2), 3));
skin_mask = logical(flip(flip(permute(niftiread('../Skin_mask.nii.gz'),[2 3 1]),2),3));
muscle_mask = logical(flip(flip(permute(niftiread('../Muscle_mask.nii.gz'),[2 3 1]),2),3));

% Calculate porosity map
skull_porosity = zeros(Nx,Ny,Nz);
skull_porosity(skull_mask) = Image(skull_mask);
skull_porosity = 1 - (skull_porosity ./ max(skull_porosity(:)));

% --------------------
% MEDIUM
% --------------------

% assign medium properties
medium.sound_speed = c_water + skin_mask*(c_skin-c_water) + muscle_mask*(c_muscle-c_water) + skull_mask.*skull_porosity*(c_skull_min-c_water) + skull_mask.*(1-skull_porosity)*(c_skull_max-c_water);
medium.density = rho_water + skin_mask*(rho_skin-rho_water) + muscle_mask*(rho_muscle-rho_water) + skull_mask.*skull_porosity*(rho_skull_min-rho_water) + skull_mask.*(1-skull_porosity)*(rho_skull_max-rho_water);
%alpha_coeff_true = alpha_coeff_water + skull_mask*(alpha_coeff_skull-alpha_coeff_water) + skin_mask*(alpha_coeff_skin-alpha_coeff_water) + muscle_mask*(alpha_coeff_muscle-alpha_coeff_water);
alpha_coeff_true =  alpha_coeff_water ...
                  + skull_mask.*((skull_porosity).^0.5).*(alpha_coeff_skull-alpha_coeff_water) ...
                  + skin_mask*(alpha_coeff_skin-alpha_coeff_water) ...
                  + muscle_mask*(alpha_coeff_muscle-alpha_coeff_water);
alpha_power_true = alpha_power_skull .* ones(Nx, Ny, Nz);

% account for actual absorption behaviour in k-Wave, which varies when high
% absorption is used (see https://doi.org/10.1121/1.4894790)
medium.alpha_coeff = fitPowerLawParamsMulti(alpha_coeff_true, alpha_power_true, medium.sound_speed, source_f0, alpha_power_skull);
%medium.alpha_coeff = alpha_coeff_true;
medium.alpha_power = alpha_power_skull;

% --------------------
% GRID
% --------------------

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);
% define axis vectors for plotting

% create the time array using an integer number of points per period
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(min(medium.sound_speed(:)) / (dx * source_f0))]);
disp(['CFL = ' num2str(max(medium.sound_speed(:)) * dt / dx)]);

% --------------------
% SOURCE
% --------------------
disp('Setting up the source ...');

% create time varying source
source_sig = createCWSignals(kgrid.t_array, source_f0, source_amp, 0);

% set bowl position and orientation
bowl_pos = [kgrid.x_vec(Bowl.x), kgrid.y_vec(Bowl.y), kgrid.z_vec(Bowl.z)];
focus_pos = [kgrid.x_vec(Focus.x), kgrid.y_vec(Focus.y), kgrid.z_vec(Focus.z)];

% create empty kWaveArray
karray = kWaveArray('SinglePrecision', true);

disp('      Adding bowl element ...');
% add bowl shaped element
karray.addBowlElement(bowl_pos, source_roc, source_diameter, focus_pos);

disp('      Creating binary mask ...');
% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

disp('      Assigning source signal ...');
% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

% --------------------
% SENSOR
% --------------------
disp('Setting up the sensor ...');

% set sensor mask to record central plane, not including the source point
sensor.mask = ones(Nx, Ny, Nz);
%sensor.mask(:, round(Ny/3):end, :) = 1;

% record the pressure
sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - round(record_periods/(source_f0*dt)) + 1;

% --------------------
% SIMULATION
% --------------------
disp('Setting up the simulation settings ...');

% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DataCast', 'single', ...
    'Smooth', [true, false, false], ...
    'SaveToDisk', 'input.h5', ...
    'DisplayMask', 'off'};

%'HDFCompressionLevel', 0, ...
%'StreamToDisk', 200, ...
%'CartInterp', 'linear', ...
%'Smooth', [true, false, false], ...
% 'DataPath' [pwd,'/Temp'], ...

grid_weights = karray.getArrayGridWeights(kgrid);
vtkwrite('Transducer_PosN.vtk', 'structured_points', 'Transducer_PosN', grid_weights);
vtkwrite('MRI_CT.vtk', 'structured_points', 'MRI_CT', Image);