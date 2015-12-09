% Author: Ke Yijie @20151209
% The program is used to solve the spatial mechanism for vectoring thrust
% in pitch channel of rotor frame
clear
clc
    
%input lengths: unit: mm
d1=20; % servo arm 
R=34.1; % gimbal radius
d3=14.25; % the y direction distance between the servo arm center and gimbal center (default: negative)
d4=20; % the x direction distance between the servo arm center and gimbal center (default: negative
l2=77; % linkage length: BC

% initial position of servo
theta1_init=95; 
servo_range=[-45,45];
stepNo=10;

% direction of point B
fai=50;  
fai=deg2rad(fai); % convert to rad

%% d5: the z direction distance between the servo arm center and gimbal
% center (default: positive), this is a redundant parameter and can be
% computed from input parameters
syms d5_sym
d5_solver=solve(R^2-2*R*d4*sin(fai)+d4^2+d3^2+d1^2-2*R*d3*cos(fai)+(2*d1*d3-2*R*d1*cos(fai))*sin(deg2rad(theta1_init))+d5_sym^2+2*d1*cos(deg2rad(theta1_init))*d5_sym-l2^2);
disp('The z directional distance between servo center and gimbal center is:')
d5=vpa(double(d5_solver(1)),8)



%% servo input and coordinates
servoinput=linspace(-1,1,stepNo);
theta1_array_deg=linspace(theta1_init-servo_range(1),theta1_init-servo_range(2),stepNo); % theta1 is ordered from big to small coresponding to servo input from -1 to 1
theta1_array=deg2rad(theta1_array_deg);

O=[0,0,0];
A=[0,-R*cos(fai),0];
D=[-d4,-d3,d5];

 %syms theta2
syms ctheta2 % cos(theta2)
% AB=[-R*sin(fai)*cos(theta2),0,-R*sin(fai)*sin(theta2)];
% 
% syms theta1
%    DC=[0,-d1*sin(theta1),-d1*cos(theta1)];
%    BC=D+DC-(A+AB)

%% 1. solve BC length equation for theta2 array
theta2_array=zeros(1, size(theta1_array,2));
theta2_array_deg=zeros(1, size(theta1_array,2));

for i=1:size(theta1_array,2)
  ctheta2_solver=vpa(solve(d4^2-2*R*d4*sin(fai)*ctheta2+R^2-2*R*d3*cos(fai)-2*R*d1*cos(fai)*sin(theta1_array(i))+2*d1*d3*sin(theta1_array(i))+d1^2+d3^2+(2*R*d5*sin(fai)-2*R*d1*sin(fai)*cos(theta1_array(i)))*sqrt(1-ctheta2^2)-2*d5*d1*cos(theta1_array(i))+d5^2-l2^2),8)
  ctheta2_pos=max(double(ctheta2_solver)); % select the positive value
  
      if(theta1_array_deg(i)<=theta1_init) % theta2>0 when the servo rotates CW and the theta1 is smaller than theta1_init
        theta2_array(i)=acos(ctheta2_pos);
        theta2_array_deg(i)=rad2deg(theta2_array(i));
      else
         theta2_array(i)=-acos(ctheta2_pos);
         theta2_array_deg(i)=rad2deg(theta2_array(i));
      end
  
end

%% 2. Polyfit the vectoring thrust pitch angle to servo input
p=polyfit(servoinput,theta2_array_deg,1);
pitch_poly=polyval(p,servoinput);

plot(servoinput,theta2_array_deg,'*',servoinput,pitch_poly,'-')
xlabel('servo input [-1,1]')
ylabel('vectoring pitch angle (deg)')
legend('computational result','polufit result')



