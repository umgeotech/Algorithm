clc;clear
% this grogram is used to calculate the excess pore water and air pressures
% during the consolidation process of unsaturated soils. This program is
% based on the governing equations proposed by Fredlund and Hasan (1979) for
% the consolidation of unsaturated soils. The readers can also refer the
% research article directly related to this grogram.
% References:
% WH Zhou, LS Zhao, XB Li. (2014) A simple analytical solution to one?dimensional 
% consolidation for unsaturated soils  International Journal for Numerical and Analytical 
% Methods in Geomechanics, 38(8), 794-810.
% Fredlund DG, Hasan JU. One-dimensional consolidation of layered systems. 
% Canadian Geotechnical Journal 1979; 17:521Ц531.
%% physical parameters
n0   = 0.5; % the initial porosity during the consolidation
Sro = 0.8;  % the degree of saturation during the consolidation
kw   = 10^(-10); % (unit:m/s)the coefficients of permeability for the water phase
ka  = 10^(-9);  % (unit:m/s)the coefficients of permeability for the air phase
ms1k = -2.5*10^(-4); % (unit:kPa-1)(ms1k=ma1k+mw1k)the volume change coefficient of soil element with respect to the net stress
ms2  = 0.4*ms1k;     % (unit:kPa-1)(ms2=ma2+mw2)the volume change coefficient of soil element with respect to the matric suction
mw1k = 0.2*ms1k;     % (unit:kPa-1)the volume change coefficient of the water phases with respect to the net stress
mw2  = 4*mw1k;       % (unit:kPa-1)the volume change coefficient of the water phases with respect to the matric suction
ma1k = ms1k - mw1k;  % (unit:kPa-1)the volume change coefficient of the air phases with respect to the net stress
ma2  = ms2  - mw2;   % (unit:kPa-1)the volume change coefficient of the air phases with respect to the matric suction
uatm = 101.3; % (unit:kPa)the atmospheric pressure
ua0  = 20; %(unit:kPa)the excess pore air pressure, it is assumed as 20 kPa at the first stage to calculate u0ahe 
u0ahe = ua0 + uatm; 
R    = 8.314; % (unit:J.mol?1.K?1) the universal air constant 
TC = 300; % (tбу+ 273)the absolute temperature (K)
molecu_mass = 29; % (unit:g.mol-1)the molecular mass of air phase 
g = 9.81; % (unit:m.seconds?2)the gravitational constant
PI = pi; % 
pw   = 1000; % (unit:kg/m3)density of water    
bata = 1;  % varied with the ratio%
Cw  = -1*(mw1k -mw2)/mw2; % pay attenttion the negative sign?    
Cwv = ms1k/mw2; 
Ca  = -1*ma2/(ma1k-ma2-uatm*n0*(1-Sro)/(u0ahe)^2);
%Cav = bata*179.3274591;
Cav = ((ka/kw)*pw*ms1k)*(R*TC/molecu_mass)/(ma2*u0ahe*(ma1k/ma2-1)-n0*(1-Sro));
%% initial distribution of excess pore pressures (water and air)
Qult = 100; % (unit:kPa)extral loading
ua0  = 20; uw0 = 40; % (unit:kPa)excess pore pressures induced by extral loading 
% ua0 is air phase uw0 is water phase
ua00 = ua0/Qult;     uw00 = uw0/Qult; % normalized excess pore pressures
%% boundary condition
% when theta is 0, it means double drainage while when theta is 1, it means
% single drainge condition.
theta = input('Please input boundary condition: Double drainge (0),single drainge (1)\n');
                
%% calculation process
Wa = Cw*Cav/(1-Ca*Cw);         Ww = Cwv/(1-Ca*Cw);
Aa = Cav/(1-Ca*Cw);            Aw = Ca*Cwv/(1-Ca*Cw);

Q1 = 0.5*(Aa+Ww+((Aa-Ww)^2+4*Aw*Wa)^0.5);        Q2 = 0.5*(Aa+Ww-((Aa-Ww)^2+4*Aw*Wa)^0.5);
Q12 = Wa/(Q2-Aa);                                Q21 = Aw/(Q1-Ww);

u10 = ua00 + Q21*uw00;
u20 = uw00 + Q12*ua00;
z  = 0.8;
H  = 1.0;

nn_point = 200;
mm = 0;
for TT = 0:1:nn_point
    T  = 10^(-8-(-9)/nn_point*TT);
    t(TT+1) = 10^(-8-(-9/nn_point)*TT);
    mm = mm+1;
       u1(mm) = 0;
       u2(mm) = 0;
    for m = 0:1:100
        if theta == 0
            M = (2*m+1)*pi;
            mid_para = 4;
        else
            M = (2*m+1)*pi/2;
            mid_para = 2;
        end
        u1(mm) = u1(mm) + mid_para*u10/M*sin(M*z/H)*exp(-M^2*T*Q1);
        u2(mm) = u2(mm) + mid_para*u20/M*sin(M*z/H)*exp(-M^2*T*Q2);
        ua(mm) = (Q21*u2(mm)-u1(mm))/(Q12*Q21-1);
        uw(mm) = (Q12*u1(mm)-u2(mm))/(Q12*Q21-1);
    end
end

%% results and plots
figure(1)
semilogx(t,ua,'-');
hold on;
figure(2)
semilogx(t,uw,'-');
hold on;

