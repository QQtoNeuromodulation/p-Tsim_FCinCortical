% Author: Soroosh Sanatkhani
% Columbia University
% Created: February 2, 2023
% Last Modified: February 6, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

load('preThermal.mat');
Water_temperature = 22; % Remove next time since it is included in part2
source_thermal.T0 = tissue.*(37-Water_temperature) + Water_temperature; % Remove next time since it is included in part2

% create kWaveDiffusion object
source_thermal.Q = Q;
kdiff = kWaveDiffusion(kgrid_thermal, medium_thermal, source_thermal, [],'DisplayUpdates', false, 'PlotSim', false);
T_t = zeros(600,length(Q(:,1,1)),length(Q(1,:,1)),length(Q(1,1,:)));

%  Pre-Cooling
kdiff.Q = 0;
dt = 5000e-6;
for i = 1:120
    disp(i);
    kdiff.takeTimeStep(round(0.5 / dt), dt);
    T_t(i,:,:,:) = kdiff.T;
end
T_min = kdiff.T;

% Sonicating
for i = 121:360
    disp(i);
    
    kdiff.Q = Q;

    dt = 100e-6;
    % take time steps
    kdiff.takeTimeStep(round(on_time / dt), dt);    
    % turn off heat source and take time steps
    kdiff.Q = 0;

    dt = 100e-6;
    kdiff.takeTimeStep(round(10e-3 / dt), dt);
    dt = 5000e-6;
    kdiff.takeTimeStep(round(off_time-20e-3 / dt), dt);
    dt = 100e-6;
    kdiff.takeTimeStep(round(10e-3 / dt), dt);
    % store the current temperature field
    
    kdiff.Q = Q;

    T_t(i,:,:,:) = kdiff.T;
end
T_max = kdiff.T;

%  Cooling
kdiff.Q = 0;
dt = 5000e-6;
for i = 361:600
    disp(i);
    kdiff.takeTimeStep(round(0.5 / dt), dt);
    T_t(i,:,:,:) = kdiff.T;
end
T_cooled = kdiff.T;

% Write VTK
vtkwrite('T_min_PosN.vtk', 'structured_points', 'T_min_PosN', T_min);
vtkwrite('T_max_PosN.vtk', 'structured_points', 'T_max_PosN', T_max);
vtkwrite('T_cooled_PosN.vtk', 'structured_points', 'T_cooled_PosN', T_cooled);
vtkwrite('Q_PosN.vtk', 'structured_points', 'Q_PosN', Q);
%
T_min_PosN = T_min;
T_max_PosN = T_max;
T_cooled_PosN = T_cooled;
Q_PosN = Q;
T_t_PosN = T_t;
% Saving
save('PosN_thermal.mat','-v7.3','T_max_PosN','T_cooled_PosN', 'T_min_PosN', 'Q_PosN','T_t_PosN');
CheckPoint = fopen('CheckPoint_Part3.txt','w');
fprintf(CheckPoint,"%s","Part 3: Successfull");
fclose(CheckPoint);
quit;