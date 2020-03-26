function rotmat=xrot(theta)
%% theta is in rad
%% left handed rotations
theta_rad=theta;
rotmat=[1 0 0; 0 cos(theta_rad) sin(theta_rad); 0 -sin(theta_rad) cos(theta_rad)];
