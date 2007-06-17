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
edge_angle = linspace(0,2*pi,n)';

% normal vectors
n1 = [cos(n1_angle),sin(n1_angle),zeros(n,1)];
n2 = [cos(n2_angle),sin(n2_angle),zeros(n,1)];

% cosine of normal angle
normal_angle_cosine = dot(n1,n2,2)./sqrt(dot(n1,n1,2))./sqrt(dot(n2,n2,2));

% normal angle
normal_angle = acos(normal_angle_cosine);

% reference vector
refvec = [zeros(n,1),zeros(n,1),-1*ones(n,1)];

% cross product of normal vectors
crossvec = cross(n1,n2,2);

% dot product of cross vector and reference vector
d = dot(crossvec,refvec,2);

% copute edge angle from normal angle and d
compute_edge_angle = zeros(n,1);
for i=1:n,
	if d(i)>0
		computed_edge_angle(i)=pi-normal_angle(i);
	else
		computed_edge_angle(i)=pi+normal_angle(i);
	end
end


% plots
figure(1);
plot(v,n1_angle*180/pi,'b',v,n2_angle*180/pi,'r',v,edge_angle*180/pi,'g');
hold on;
plot(v,normal_angle*180/pi,'c',v,computed_edge_angle*180/pi,'k');
hold off;
legend('n1 angle','n2 angle','edge angle','normal angle','compute edge angle');

figure(2);
plot(v,(pi-normal_angle)*180/pi,'b',v,(pi+normal_angle)*180/pi,'r');
legend('pi-normal angle','pi+normal angle');

figure(3);
plot(v,d);
legend('d');

