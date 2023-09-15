% Likelihood function
% 3 unknown parameters
% The measurement data (TNEC_deflection) can be found in the paper of Ou et al

function log_like =TNEC_post_3(x,exponential)
strut1=0*2+1;
strut2=35*2+1;
strut3=71*2+1; 
strut4=103*2+1;
strut5=137*2+1;
strut6=165*2+1;
impact_depth=4;
m=350;
depth_indice=0:0.1:34.9;		
b=depth_indice';
TNEC_deflection_strut1=20/1000;
TNEC_deflection_strut2=29.47/1000;
TNEC_deflection_strut3=40.08/1000;
TNEC_deflection_strut4=56.74/1000;
TNEC_deflection_strut5=75.48/1000;
TNEC_deflection_strut6=98.58/1000;

EI=1275750;% stiffness of the wall stiffness
s1=152913;%stiffness of the concrete slab stiffness
cf=x(2); % soil spring constant
a=x(3);%slope of soil pressure above surface
L=35;
data_size =35;
%% stage 3
nodeCoordinates=linspace(0,L,m+1)';%L is the length of beam, constructing a column matrix
numberNodes=m+1;
GDof=2*numberNodes;%Number for degrees of freedom
U=zeros(GDof,1); %Solution for the equation
force=zeros(GDof,1); %force matrix
stiffness=zeros(GDof); %stiffness matrix
for i=1:m;
elementNodes(i,1)=i;
elementNodes(i,2)=i+1;%elementNodes matrix
end
for j=1:m;
indice=elementNodes(j,:);%indice matrix using for constructing stiffness matrix
elementIndice=[2*(indice(1)-1)+1 2*(indice(2)-1) 2*(indice(2)-1)+1 2*(indice(2)-1)+2];
LElem=nodeCoordinates(indice(2))-nodeCoordinates(indice(1));%Calculating distance of two adjacent element
kk=2*EI/(LElem)^3*[6 -3*LElem -6 -3*LElem;-3*LElem 2*LElem^2 3*LElem LElem^2; -6 3*LElem 6 3*LElem; -3*LElem LElem^2 3*LElem 2*LElem^2];%Constructing the stiffness matrix of EI
stiffness(elementIndice,elementIndice)=stiffness(elementIndice,elementIndice)+kk;%assembly of stiffness matrix according to indice
end
excavation_depth=8.6;
q=a*excavation_depth;
for j=1:86;
indice=elementNodes(j,:);%indice matrix using for constructing stiffness matrix
elementIndice=[2*(indice(1)-1)+1 2*(indice(2)-1) 2*(indice(2)-1)+1 2*(indice(2)-1)+2];
LElem=nodeCoordinates(indice(2))-nodeCoordinates(indice(1));%Calculating distance of two adjacent element
ff=[(a*LElem*(10*b(j,:)+3*LElem))/20 -(a*LElem^2*(5*b(j,:)+2*LElem))/60  (a*LElem*(10*b(j,:)+7*LElem))/20  (a*LElem^2*(5*b(j,:)+3*LElem))/60]';
force(elementIndice)=force(elementIndice)+ff; 
end
for k=87:m;
indice=elementNodes(k,:);%indice matrix using for constructing stiffness matrix
elementIndice=[2*(indice(1)-1)+1 2*(indice(2)-1) 2*(indice(2)-1)+1 2*(indice(2)-1)+2];
LElem=nodeCoordinates(indice(2))-nodeCoordinates(indice(1));%Calculating distance of two adjacent element
ff=q*LElem/12*[6 -LElem 6 LElem]';%Constructing the force matrix
force(elementIndice)=force(elementIndice)+ff; %assembly of uniformly distributed load matrix
end
for k=87:m;
indice=elementNodes(k,:);%indice matrix using for constructing stiffness matrix
elementIndice=[2*(indice(1)-1)+1 2*(indice(2)-1) 2*(indice(2)-1)+1 2*(indice(2)-1)+2];
LElem=nodeCoordinates(indice(2))-nodeCoordinates(indice(1));%Calculating distance of two adjacent element
c= b(k-86,1);
switch exponential
case 1  %n=0
s11 =(13*LElem*cf)/35;
s22 =(LElem^3*cf)/105;
s33 =(13*LElem*cf)/35;
s44 =(LElem^3*cf)/105;
s12 =-(11*LElem^2*cf)/210;
s13 =(9*LElem*cf)/70;
s14 =(13*LElem^2*cf)/420;
s23 =-(13*LElem^2*cf)/420;
s24 =-(LElem^3*cf)/140;
s34 =(11*LElem^2*cf)/210;
case 2 %n=0.5
s11 =(cf*(LElem + c)^(3/2)*((7424*LElem^6)/45045 + (1024*LElem^5*c)/3003 - (512*LElem^4*c^2)/5005 - (23552*LElem^3*c^3)/45045 - (256*LElem^2*c^4)/15015 + (2048*LElem*c^5)/5005 + (8192*c^6)/45045))/LElem^6 - (c^(3/2)*cf*((2*LElem^6)/3 - (32*LElem^4*c^2)/35 - (128*LElem^3*c^3)/315 + (256*LElem^2*c^4)/385 + (2048*LElem*c^5)/3003 + (8192*c^6)/45045))/LElem^6;
s22 =(cf*((LElem + c)^(7/2)*((2*LElem^4)/7 + (24*LElem^3*c)/7 + (72*LElem^2*c^2)/7 + (80*LElem*c^3)/7 + (30*c^4)/7) + c^(5/2)*((4*LElem^4*c)/5 + (24*LElem^3*c^2)/5 + (48*LElem^2*c^3)/5 + 8*LElem*c^4 + (12*c^5)/5) + (LElem + c)^(3/2)*((2*LElem^4*c^2)/3 + (8*LElem^3*c^3)/3 + 4*LElem^2*c^4 + (8*LElem*c^5)/3 + (2*c^6)/3) - (LElem + c)^(9/2)*((8*LElem^3)/9 + (16*LElem^2*c)/3 + (80*LElem*c^2)/9 + (40*c^3)/9) - c^(7/2)*((2*LElem^4)/7 + (24*LElem^3*c)/7 + (72*LElem^2*c^2)/7 + (80*LElem*c^3)/7 + (30*c^4)/7) + (LElem + c)^(11/2)*((12*LElem^2)/11 + (40*LElem*c)/11 + (30*c^2)/11) + (2*(LElem + c)^(15/2))/15 - c^(3/2)*((2*LElem^4*c^2)/3 + (8*LElem^3*c^3)/3 + 4*LElem^2*c^4 + (8*LElem*c^5)/3 + (2*c^6)/3) + c^(9/2)*((8*LElem^3)/9 + (16*LElem^2*c)/3 + (80*LElem*c^2)/9 + (40*c^3)/9) - (LElem + c)^(13/2)*((8*LElem)/13 + (12*c)/13) - c^(11/2)*((12*LElem^2)/11 + (40*LElem*c)/11 + (30*c^2)/11) - (2*c^(15/2))/15 - (LElem + c)^(5/2)*((4*LElem^4*c)/5 + (24*LElem^3*c^2)/5 + (48*LElem^2*c^3)/5 + 8*LElem*c^4 + (12*c^5)/5) + c^(13/2)*((8*LElem)/13 + (12*c)/13)))/LElem^4;
s33 =(694*cf*(LElem + c)^(3/2))/2145 + ((2*cf*(4096*c^6*(LElem + c)^(3/2) - 4096*c^(15/2)))/45045 + (2*LElem*cf*(9216*c^5*(LElem + c)^(3/2) - 15360*c^(13/2)))/45045 - (2*LElem^2*cf*(384*c^4*(LElem + c)^(3/2) + 14976*c^(11/2)))/45045 - (192*LElem^5*c*cf*(LElem + c)^(3/2))/715 - (5248*LElem^3*c^3*cf*(LElem + c)^(3/2))/45045 + (608*LElem^4*c^2*cf*(LElem + c)^(3/2))/3003)/LElem^6;
s44 =-(cf*((1024*LElem*c^(13/2))/9009 - (16*LElem^6*(LElem + c)^(3/2))/2145 - (2048*c^6*(LElem + c)^(3/2))/45045 + (2048*c^(15/2))/45045 + (256*LElem^2*c^(11/2))/3465 + (512*LElem^2*c^4*(LElem + c)^(3/2))/45045 - (128*LElem^3*c^3*(LElem + c)^(3/2))/45045 - (16*LElem^4*c^2*(LElem + c)^(3/2))/9009 - (2048*LElem*c^5*(LElem + c)^(3/2))/45045 + (32*LElem^5*c*(LElem + c)^(3/2))/6435))/LElem^4;
s12 =-(12012*LElem^5*c^(5/2)*cf - 4096*c^(15/2)*cf + 13728*LElem^4*c^(7/2)*cf - 9152*LElem^3*c^(9/2)*cf - 26624*LElem^2*c^(11/2)*cf - 17920*LElem*c^(13/2)*cf + 1280*LElem^6*cf*(LElem + c)^(3/2) + 4096*c^6*cf*(LElem + c)^(3/2) + 11776*LElem*c^5*cf*(LElem + c)^(3/2) + 512*LElem^5*c*cf*(LElem + c)^(3/2) + 7424*LElem^2*c^4*cf*(LElem + c)^(3/2) - 6144*LElem^3*c^3*cf*(LElem + c)^(3/2) - 6656*LElem^4*c^2*cf*(LElem + c)^(3/2))/(45045*LElem^5);
s13 =(8192*c^(15/2)*cf - 20592*LElem^4*c^(7/2)*cf - 9152*LElem^3*c^(9/2)*cf + 29952*LElem^2*c^(11/2)*cf + 30720*LElem*c^(13/2)*cf + 4016*LElem^7*cf*(LElem + c)^(1/2) - 8192*c^7*cf*(LElem + c)^(1/2) - 26624*LElem*c^6*cf*(LElem + c)^(1/2) + 2384*LElem^6*c*cf*(LElem + c)^(1/2) - 17664*LElem^2*c^5*cf*(LElem + c)^(1/2) + 15168*LElem^3*c^4*cf*(LElem + c)^(1/2) + 12144*LElem^4*c^3*cf*(LElem + c)^(1/2) - 3888*LElem^5*c^2*cf*(LElem + c)^(1/2))/(45045*LElem^6);
s14 =(4096*c^(15/2)*cf - 6864*LElem^4*c^(7/2)*cf - 4576*LElem^3*c^(9/2)*cf + 9984*LElem^2*c^(11/2)*cf + 12800*LElem*c^(13/2)*cf + 928*LElem^7*cf*(LElem + c)^(1/2) - 4096*c^7*cf*(LElem + c)^(1/2) - 10752*LElem*c^6*cf*(LElem + c)^(1/2) + 640*LElem^6*c*cf*(LElem + c)^(1/2) - 5120*LElem^2*c^5*cf*(LElem + c)^(1/2) + 6048*LElem^3*c^4*cf*(LElem + c)^(1/2) + 3712*LElem^4*c^3*cf*(LElem + c)^(1/2) - 1088*LElem^5*c^2*cf*(LElem + c)^(1/2))/(45045*LElem^5);
s23 =-(4096*c^(15/2)*cf + 13728*LElem^3*c^(9/2)*cf + 26624*LElem^2*c^(11/2)*cf + 17920*LElem*c^(13/2)*cf + 1008*LElem^6*cf*(LElem + c)^(3/2) - 4096*c^6*cf*(LElem + c)^(3/2) - 11776*LElem*c^5*cf*(LElem + c)^(3/2) - 512*LElem^5*c*cf*(LElem + c)^(3/2) - 7424*LElem^2*c^4*cf*(LElem + c)^(3/2) + 1568*LElem^3*c^3*cf*(LElem + c)^(3/2) - 208*LElem^4*c^2*cf*(LElem + c)^(3/2))/(45045*LElem^5);
s24 =-(2048*c^(15/2)*cf + 4576*LElem^3*c^(9/2)*cf + 9984*LElem^2*c^(11/2)*cf + 7680*LElem*c^(13/2)*cf + 224*LElem^6*cf*(LElem + c)^(3/2) - 2048*c^6*cf*(LElem + c)^(3/2) - 4608*LElem*c^5*cf*(LElem + c)^(3/2) - 96*LElem^5*c*cf*(LElem + c)^(3/2) - 2304*LElem^2*c^4*cf*(LElem + c)^(3/2) + 480*LElem^3*c^3*cf*(LElem + c)^(3/2) - 96*LElem^4*c^2*cf*(LElem + c)^(3/2))/(45045*LElem^4);
s34 =(2*cf*((2*(LElem + c)^(1/2)*(3003*LElem^7 + 231*LElem^6*c - 252*LElem^5*c^2 + 280*LElem^4*c^3 - 320*LElem^3*c^4 + 384*LElem^2*c^5 - 512*LElem*c^6 + 1024*c^7))/45045 - (2048*c^(15/2))/45045))/LElem^5 - (5*cf*((512*c^(13/2))/9009 + (2*(LElem + c)^(3/2)*(693*LElem^5 - 630*LElem^4*c + 560*LElem^3*c^2 - 480*LElem^2*c^3 + 384*LElem*c^4 - 256*c^5))/9009))/LElem^4 + (3*cf*((2*(LElem + c)^(3/2)*(315*LElem^4 - 280*LElem^3*c + 240*LElem^2*c^2 - 192*LElem*c^3 + 128*c^4))/3465 - (256*c^(11/2))/3465))/LElem^3;
case 3 %n=1
s11 =(LElem*cf*(3*LElem+13*c))/35;
s22 =(LElem^3*cf*(3*LElem+8*c))/840;
s33 =(LElem*cf*(10*LElem+13*c))/35;
s44 =(LElem^3*cf*(5*LElem+8*c))/840;
s12 =-(LElem^2*cf*(7*LElem+22*c))/420;
s13 =(9*LElem*cf*(LElem+2*c))/140;
s14 =(LElem^2*cf*(6*LElem+13*c))/420;
s23 =-(LElem^2*cf*(7*LElem+13*c))/420;
s24 =-(LElem^3*cf*(LElem+2*c))/280;
s34 =(LElem^2*cf*(15*LElem+22*c))/420;
case 4 %n=2
s11 =(LElem*cf*(19*LElem^2+108*LElem*c+234*c^2))/630;
s22 =(LElem^3*cf*(2*LElem^2+9*LElem*c+12*c^2))/1260;
s33 =(LElem*cf*(145*LElem^2+360*LElem*c+234*c^2))/630;
s44 =(LElem^3*cf*(5*LElem^2+15*LElem*c+12*c^2))/1260;
s12 =-(LElem^2*cf*(17*LElem^2+84*LElem*c+132*c^2))/2520;
s13 =(LElem*cf*(23*LElem^2+81*LElem*c+81*c^2))/630;
s14 =(LElem^2*cf*(19*LElem^2+72*LElem*c+78*c^2))/2520;
s23 =-(LElem^2*cf*(25*LElem^2+84*LElem*c+78*c^2))/2520;
s24 =-(LElem^3*cf*(5*LElem^2+18*LElem*c+18*c^2))/2520;
s34 =(LElem^2*cf*(65*LElem^2+180*LElem*c+132*c^2))/2520;
case 5
if c<impact_depth;
%n=1
s11 =(LElem*cf*(3*LElem+13*c))/35;
s22 =(LElem^3*cf*(3*LElem+8*c))/840;
s33 =(LElem*cf*(10*LElem+13*c))/35;
s44 =(LElem^3*cf*(5*LElem+8*c))/840;
s12 =-(LElem^2*cf*(7*LElem+22*c))/420;
s13 =(9*LElem*cf*(LElem+2*c))/140;
s14 =(LElem^2*cf*(6*LElem+13*c))/420;
s23 =-(LElem^2*cf*(7*LElem+13*c))/420;
s24 =-(LElem^3*cf*(LElem+2*c))/280;
s34 =(LElem^2*cf*(15*LElem+22*c))/420;
else 
%n=0
cf_q=cf*impact_depth;
s11 =(13*LElem*cf_q)/35;
s22 =(LElem^3*cf_q)/105;
s33 =(13*LElem*cf_q)/35;
s44 =(LElem^3*cf_q)/105;
s12 =-(11*LElem^2*cf_q)/210;
s13 =(9*LElem*cf_q)/70;
s14 =(13*LElem^2*cf_q)/420;
s23 =-(13*LElem^2*cf_q)/420;
s24 =-(LElem^3*cf_q)/140;
s34 =(11*LElem^2*cf_q)/210;
end
end
s21=s12;
s31=s13;
s41=s14;
s32=s23;
s42=s24;
s43=s34;
cc=[s11,s12,s13,s14; s21,s22,s23,s24; s31,s32,s33,s34; s41,s42,s43,s44];
stiffness(elementIndice,elementIndice)=stiffness(elementIndice,elementIndice)+cc;%assembly of stiffness matrix according to indice
end
stiffness(strut1, strut1)= stiffness(strut1, strut1)+s1;
force(strut1)=force(strut1)+s1*TNEC_deflection_strut1;
stiffness(strut2, strut2)= stiffness(strut2, strut2)+s1;
force(strut2)=force(strut2)+s1*TNEC_deflection_strut2;
U=stiffness\force;
deflection=1000.*U(21:20:701);
load('TNEC_deflection.mat') 
measured_deflection=TNEC_deflection(:,3); 
log_like3=-0.5*data_size*log(2*pi)-data_size*log(x(1))-0.5*(deflection-measured_deflection)'*(deflection-measured_deflection)/x(1)^2;
log_like=log_like3;