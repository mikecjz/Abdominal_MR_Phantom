function rotmat=yrot(theta)
%% theta is in rad
%% left handed rotations
theta_rad=theta;
rotmat=[cos(theta_rad) 0 -sin(theta_rad);0 1 0; sin(theta_rad) 0 cos(theta_rad)];
