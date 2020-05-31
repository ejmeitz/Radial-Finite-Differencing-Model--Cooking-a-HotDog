clc;
%properties
C1 = 1.0712;
k = 0.52;     %thermal conductivity
roe = 880;   %density of hot dog
c = 3350;  %specific heat
zeta = 0.7465;
alpha = k/(roe*c);    %thermal diffusivity
r0 = 0.0127; %[m]

%arrays
%note the radial array is *10000 so that the indicies can still be whole
%numbers
theta_star = zeros(128,1000);   %[r*10000 , seconds]
time_array = zeros(1,1000);
radius_array = zeros(1,128);
T= zeros(128,1000);  %[theta rearranged for just T]
Ti = 10+273;      %hot dog initial temp
T_inf = 250+273;    %air temp around hot dog
 
%h_total = h_rad + h_conv
h = 7.6+5.55;
Bi = h*r0/k;

%loop through time and radial position, adjust time and radial positions
%down 1 to start at 0 within loop
for t = 1:1:1000       
    time = t-1;
    time_array(t) = time;
    for r = 1:1:128
        radius = (r-1)/10000;  %start r at 1 and put in meters
        radius_array(r) = radius;
        %calculate theta start at every radial position and time 
        theta_star(r,t) = C1*exp((-1*(zeta^2)*((alpha * time))/(r0^2)))*besselj(0,(zeta*radius/r0));
        %solve for T from theta star
        T(r,t) = (theta_star(r,t) * (Ti -T_inf)) +T_inf; 
    end
end

%create arrays to hold times
time_0 = zeros(1,128);
time_30 = zeros(1,128);
time_120 = zeros(1,128);
time_570 = zeros(1,128);

%fill in arrays
for r = 1:1:128
    time_0(1,r) = T(r,1);
    time_30(1,r) = T(r,31);
    time_120(1,r) = T(r,121);
    time_570(1,r) = T(r,571);
    
end

%plot data
figure(1); hold on;
%t0 = plot(radius_array,time_0); L1 = "t = 0s";
t30 = plot(radius_array,time_30); L2 = "t = 30s";
t120 = plot(radius_array,time_120); L3 = "t = 120s";
t306 = plot(radius_array,time_570); L5 = "t = 570s";
legend([t30,t120,t306],[L2,L3,L5]);
xlabel("Radius [m]")
ylabel("Temperature [K]")
%title("Radial Temperature at Various Times");


centerline = zeros(1,1000);
surface = zeros(1,1000);
for t = 1:1:1000
   centerline(1,t) = T(1,t);
   surface(1,t) = T(128,t);
end


figure(2); hold on;
c = plot(time_array,centerline); L6 = "Centerline Temperature";
s = plot(time_array,surface); L7 = "Surface Temperature";
legend([c,s],[L6,L7]);
xlabel("Time [sec]");
ylabel("Temperature [K]");
%title("Temperature of the hotdog at the centerline and surface with respect to time");


