clear all;
close all;
clc;

% number of samples
n = 100;
v = [1:n];

% angle of n1, radians
n1_angle = pi/2*ones(n,1);

% angle of n2, radians
n2_angle = linspace(-pi/2,3*pi/2,n)';

% edge angle, radians
theta = linspace(0,2*pi,n)';

% normal vectors
n1 = [cos(n1_angle),sin(n1_angle),zeros(n,1)];
n2 = [cos(n2_angle),sin(n2_angle),zeros(n,1)];

v1 = [1 1 1];
v2 = [1 1 10];
o1 = [3 1 3];
o2 = [2*cos(theta)+1 2*sin(theta)+1 7*ones(n,1)];

% reference vector
refvec = [repmat(v2,[n,1])-o2];

% compute edge stretch angle gamma
gamma = -dot(n1,refvec,2)./abs(dot(n1,refvec,2)).*acos(dot(n1,n2,2));
%gamma = acos(dot(n1,n2,2));

% cosine of normal angle
%normal_angle_cosine = dot(n1,n2,2)./sqrt(dot(n1,n1,2))./sqrt(dot(n2,n2,2));

% normal angle
%normal_angle = acos(normal_angle_cosine);

% cross product of normal vectors
%crossvec = cross(n1,n2,2);

% dot product of cross vector and reference vector
%d = dot(crossvec,refvec,2);

% copute edge angle from normal angle and d
%compute_edge_angle = zeros(n,1);
%for i=1:n,
%	if d(i)>0
%		computed_edge_angle(i)=pi-normal_angle(i);
%	else
%		computed_edge_angle(i)=pi+normal_angle(i);
%	end
%end


% plots
figure(1);
plot(v,n1_angle*180/pi,'b',v,n2_angle*180/pi,'r',v,theta*180/pi,'g',v,gamma*180/pi,'k');
%hold on;
%plot(v,normal_angle*180/pi,'c',v,computed_edge_angle*180/pi,'k');
%hold off;
legend('n1 angle','n2 angle','edge angle','gamma');

%figure(2);
%plot(v,(pi-normal_angle)*180/pi,'b',v,(pi+normal_angle)*180/pi,'r');
%legend('pi-normal angle','pi+normal angle');

%figure(3);
%plot(v,d);
%legend('d');

