clc; 
clear all;
close all;

%number of nodes
M = 25;
%length of heating--arbitrary just long enough to see behavior
dur = 1000;   %[s]

%properties
d = 0.0254;  %[m]  dia of hotdog
r = d/2;    %[m]  radius of hotdog
roe = 880;  %[kg/m^3]    density of hot dog
k = 0.52;  %[W/mK]       thermal conductivity of hotdog
c = 3350; %[J/kgK]       specific heat of hotdog
alpha = k/(roe*c); 

%solved for
h = 7.61;   %convective heat transfer coefficient
h_r = 5.55;  %radiative heat transfer coefficient
h = h+h_r;   %total heat transfer coefficient
delta_r = r/M;     %differential r
Bi = h*delta_r/k;         %FDE Biot number


%Stability Calculation
%The central node gives the smallest Fo which means the smallest time step
%so we use that as the Fourier number
Fo = 0.25;
time_step = ((Fo*(delta_r)^2)/alpha);

%this is the number of samples that will be taken to acheive the chosen
%duration
samples = ceil(dur/time_step);

%initial temp
T_0 = 283;  %[K]
%air temperature (infinity)
T_inf = 250+273;

%Note every 1 step in t is 1 time step not 1 second 
T = zeros(M,samples);

%instantiate initial temps
for j = 1:M+1
    T(j,1) = T_0;
end

for t = 1:samples+1
    for m = 1:M+1
        %NOTE t and m START AT 1 BECAUSE MATLAB IS STUPID, so m-1 was
        %subbed in for all m's in the equation

        % center node equation 
        if m == 1
            T(m,t+1) = 4*Fo*T(m+1,t)+(T(m,t)*(1-4*Fo));
            
        % exterior node equation
        elseif m == M+1
            %this commented out one uses V = M*dr^2 * pi * 0.5 * H
            %they're roughly the same result
            %T(m,t+1) = (Fo*(1-(1/(2*(m-1))))*(T(m-1,t)-T(m,t)))+ 2*(Bi*Fo*(T_inf-T(m,t)))+ T(m,t);
           
            %this one uses V = (pi*dr^2*(M-0.25)) *H   r_o^2 - r_i^2
            T(m,t+1) = (Fo*((m-1)-0.5)/((m-1)-0.25))*(T(m-1,t)-T(m,t)) + ((2*Bi*Fo*(m-1))/((m-1)-0.25))*(T_inf - T(m,t)) + T(m,t);
        
            %interior node equation
        else
             T(m,t+1) = Fo*(1-(1/(2*(m-1))))*(T(m-1,t)) + Fo*(1+(1/(2*(m-1))))*(T(m+1,t))+(1-(2*Fo))*(T(m,t));
        end
        
    end
end


figure(1);
hold on;
center = plot(T(1,:)); L1 = "Centerline";
surf = plot(T(M+1,:)); L2 = "Surface";
y1 = yline(68+273,'g'); L3 = "Ideal Centerline Temp";
y2 = yline(273+100,'m'); L4 = "Maximum Surface Temp";

legend([center,surf,y1,y2],[L1,L2,L3,L4]);
xlabel('Time [s]');
ylabel('Temperature [K]');

%adjust scaling on x axis
xticks = get(gca,'xtick'); 
scaling  = time_step; 
newlabels = arrayfun(@(x) sprintf('%.1f', scaling * x), xticks, 'un', 0);
set(gca,'xticklabel',newlabels);

figure(2);
hold on;
%note the positions in the temperature array are at 30/timestep since the
%temperature array is not per every second its per every time step
%t_0 = plot(1:1:M+1,T(:,1)); T1 ="t=0";
t_30 = plot(1:1:M+1,T(:,85)); T2 ="t=30";
t_90 = plot(1:1:M+1,T(:,326)); T3 ="t=120";
t_final = plot(1:1:M+1,T(:,958)); T4="Final Time = 353s";
xlabel('Radial Node #');
ylabel('Temperature [K]');
legend([t_30,t_90,t_final],[T2,T3,T4]);

