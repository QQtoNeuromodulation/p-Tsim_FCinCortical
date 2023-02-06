% Author: Soroosh Sanatkhani
% Columbia University
% Created: February 2, 2023
% Last Modified: February 6, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

load('pre.mat');
sensor_data.p = h5read('output.h5', '/p');

grid_weights = karray.getArrayGridWeights(kgrid);
clear karray;

disp('Extracting pressure ...');
pressure = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
clear kgrid;
pressure = reshape(pressure, [Nx,Ny,Nz]);

disp(['Resolution (mm) = (', num2str(resolution),' ', num2str(resolution),' ', num2str(resolution),')']);
disp(['Time Step (ns) = ', num2str(dt)]);
disp(['Simulation Length (micro second) = ', num2str(t_end*1e6)]);
disp(['PPW = ' num2str(min(medium.sound_speed(:)) / (dx * source_f0))]);
disp(['CFL = ' num2str(max(medium.sound_speed(:)) * dt / dx)]);
disp(['Maximum Pressure (kPa) = ', num2str(max(pressure(:))./1000)]);

close all;
Pressure_PosN = pressure;
sensor_data_PosN = sensor_data;
Transducer_PosN = grid_weights;
clear sensor_data grid_weights;

% Write VTK

vtkwrite('Transducer_PosN.vtk', 'structured_points', 'Transducer_PosN', Transducer_PosN);
vtkwrite('Pressure_PosN.vtk', 'structured_points', 'Pressure_PosN', Pressure_PosN);

% Extra information
vtkwrite('MRI_CT.vtk', 'structured_points', 'MRI_CT', Image);
vtkwrite('Attenuation.vtk', 'structured_points', 'Attenuation', medium.alpha_coeff);
vtkwrite('Density.vtk', 'structured_points', 'Density', medium.density);
vtkwrite('sound_speed.vtk', 'structured_points', 'sound_speed', medium.sound_speed);
vtkwrite('Skull_mask.vtk', 'structured_points', 'Skull_mask', skull_mask);

save('PosN_acoustic.mat','-v7.3','Pressure_PosN','sensor_data_PosN','Transducer_PosN','Picked_Slice');
clearvars -except kgrid_thermal medium pressure source_thermal Q on_time off_time Picked_Slice source_f0 skin_mask muscle_mask skull_mask brain_mask dx Nx Ny Nz;

% convert the absorption coefficient to nepers/m
alpha_coeff_water_absorp     = 0.0022;
alpha_coeff_skull_absorp     = 2.7;
alpha_coeff_skin_absorp      = 1.84;
alpha_coeff_muscle_absorp    = 0.62;
alpha_coeff_brain_absorp     = 0.59;
alpha_power_skull     = 1.1;
alpha_coeff_true_absorp = alpha_coeff_water_absorp ...
                        + skull_mask*(alpha_coeff_skull_absorp-alpha_coeff_water_absorp) ...
                        + skin_mask*(alpha_coeff_skin_absorp-alpha_coeff_water_absorp) ...
                        + muscle_mask*(alpha_coeff_muscle_absorp-alpha_coeff_water_absorp) ...
                        + brain_mask*(alpha_coeff_brain_absorp-alpha_coeff_water_absorp);
alpha_power_true_absorp = alpha_power_skull .* ones(Nx, Ny, Nz);
% account for actual absorption behaviour in k-Wave, which varies when high
% absorption is used (see https://doi.org/10.1121/1.4894790)
medium.alpha_coeff_absorp = fitPowerLawParamsMulti(alpha_coeff_true_absorp, alpha_power_true_absorp, medium.sound_speed, source_f0, alpha_power_skull);

alpha_np = db2neper(medium.alpha_coeff_absorp, medium.alpha_power) * ...
    (2 * pi * source_f0).^medium.alpha_power;

Q = alpha_np .* pressure.^2 ./ (medium.density .* medium.sound_speed);

% clear the input structures
%clear medium source sensor;
% set the background temperature and heating term

% set source on time and off time %%% 240 pulses; on = 10 ms; off = 490 ms;
% 2 minutes total
on_time  = 10e-3;  % [s]
off_time = 490e-3;  % [s]

% set time step size

Cp_water = 4178;
Cp_skin = 3391;
Cp_muscle = 3421;
Cp_brain = 3630;
Cp_skull = 1313;
Cp_blood = 3617;

k_water = 0.6;
k_skin = 0.37;
k_muscle = 0.49;
k_brain = 0.51;
k_skull = 0.3;

rho_skin        = 1109; % density [kg/m^3]
rho_muscle      = 1090;
rho_brain       = 1046;
rho_blood       = 1050;

w_muscle = 0.0006; %[1/s] (actually brain) (37 mL/(min.kg))
w_brain = 0.0098;
w_skin = 0.0019; %(106 mL/(min.kg)) [(106*1e-6*1050)/60]
%w_skull = 0;

% define medium properties related to diffusion
medium_thermal.density = medium.density;
medium_thermal.thermal_conductivity = k_water + skin_mask*(k_skin-k_water) + muscle_mask*(k_muscle-k_water) + brain_mask*(k_brain-k_water) + skull_mask.*(k_skull-k_water);
medium_thermal.specific_heat        = Cp_water + skin_mask*(Cp_skin-Cp_water) + muscle_mask*(Cp_muscle-Cp_water) + brain_mask*(Cp_brain-Cp_water) + skull_mask.*(Cp_skull-Cp_water);     % [J/(kg.K)]
medium_thermal.perfusion_coeff      = 0 + skin_mask.*((rho_blood.*w_skin.*Cp_blood)./(rho_skin.*Cp_skin)) ...
                                        + muscle_mask.*((rho_blood.*w_muscle.*Cp_blood)./(rho_muscle.*Cp_muscle)) ...
                                        + brain_mask.*((rho_blood.*w_brain.*Cp_blood)./(rho_brain.*Cp_brain)) ...
                                        + 0;
% kdiff.perfusion_coeff =   medium.blood_density .*
%                           medium.blood_perfusion_rate .* medium.blood_specific_heat ./
%                           (medium.density .* medium.specific_heat);

vtkwrite('thermal_conductivity.vtk', 'structured_points', 'thermal_conductivity', medium_thermal.thermal_conductivity);
vtkwrite('specific_heat.vtk', 'structured_points', 'specific_heat', medium_thermal.specific_heat);

%%Cropping for thermal
Q = Q(81:317,121:290,51:251);

medium_thermal.thermal_conductivity = medium_thermal.thermal_conductivity(81:317,121:290,51:251);
medium_thermal.specific_heat = medium_thermal.specific_heat(81:317,121:290,51:251);
medium_thermal.density = medium_thermal.density(81:317,121:290,51:251);
medium_thermal.perfusion_coeff = medium_thermal.perfusion_coeff(81:317,121:290,51:251);
tissue = skin_mask(81:317,121:290,51:251)+muscle_mask(81:317,121:290,51:251)+brain_mask(81:317,121:290,51:251)+skull_mask(81:317,121:290,51:251);

Water_temperature = 22;
medium_thermal.blood_ambient_temperature = 37;
source_thermal.T0 = tissue.*(37-Water_temperature) + Water_temperature;

% create the computational grid

[Nx, Ny, Nz] = size(Q);
kgrid_thermal = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

clearvars -except kgrid_thermal medium_thermal source_thermal Water_temperature Q on_time off_time Picked_Slice tissue;
save('preThermal.mat','-v7.3');
CheckPoint = fopen('CheckPoint_Part2.txt','w');
fprintf(CheckPoint,"%s","Part 2: Successfull");
fclose(CheckPoint);
quit;