clc;
close all;
clear all;


%properties
e_coals = 0.8; %emissivity of coals
e_HD = 0.35; %emissivity of HD
s1 = 0.2032;  %[m]
s2 = -0.2030;  %[m]
L = 0.1016;   %[m]
r_hd = 0.0127;  % [m] radius of HD
T_coals = 450+273;  %[K] temperature of coals


%calculation for view factor

F_ij = (r_hd/(s1-s2))*(atan(s1/L)-atan(s2/L));
A_i = 0.4826;   % area of coals per unit length  (s1-s2)
A_j = 0.079;  %area of hotdog per unit length    (circumference)
F_ji = (A_i*F_ij)/A_j;

T2_star = (60+273)/((e_coals*F_ji)^0.25);
h_rad = e_coals*e_HD*F_ji*(5.67*10^(-8))*(T_coals+T2_star)*(((T_coals^2)+(T2_star)^2));