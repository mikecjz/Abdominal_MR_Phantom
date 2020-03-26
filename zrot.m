function rotmat=zrot(theta)
%% theta is in rad
%% left handed rotations
theta_rad=theta;
rotmat=[cos(theta_rad) sin(theta_rad) 0; -sin(theta_rad) cos(theta_rad) 0; 0 0 1];
