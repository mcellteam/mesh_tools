clear all;
close all;
clc;

% number of samples
n = 100;
v = [1:n];

% angle of edge, radians
angle = linspace(0,2*pi,n)';

v1 = [1 1 1];
v2 = [1 1 10];
o1 = [3 1 3];
o2 = [2*cos(angle)+1 2*sin(angle)+1 7*ones(n,1)];

% find perpendicular intersection points
delv = v2-v1;
%den = dot(v1,v1)-2*dot(v1,v2)+dot(v2,v2);
den = dot(delv,delv,2);
if den~=0
    num = dot(delv,o1,2)-dot(delv,v1);
	i1 = v1+num/den*(v2-v1);
    num = dot(repmat(delv,[n,1]),o2,2)-dot(delv,v1);
    %num = dot(v11,v2)*one(v1repmat(v2,[n,1]),o2,2);
	i2 = repmat(v1,[n,1])+num./repmat(den,[n,1])*(v2-v1);
end


% compute clock hands
s1 = o1-i1;
s2 = o2-i2;

% compute angular spring stretch
gamma_first = acos(-dot(repmat(s1,[n,1]),s2,2)./sqrt(dot(repmat(s1,[n,1]),repmat(s1,[n,1]),2))./sqrt(dot(s2,s2,2)));

gamma = [gamma_first(1:50);-gamma_first(51:100)];
% cosine of normal angle
%normal_angle_cosine = dot(n1,n2,2)./sqrt(dot(n1,n1,2))./sqrt(dot(n2,n2,2));

% normal angle
%normal_angle = acos(normal_angle_cosine);

% reference vector
%refvec = [zeros(n,1),zeros(n,1),-1*ones(n,1)];

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
%

% plots
figure(1);
plot(v,angle*180/pi,'b',v,gamma*180/pi,'r',v,(gamma+angle)*180/pi,'c');
%hold on;
%plot(v,normal_angle*180/pi,'c',v,computed_edge_angle*180/pi,'k');
%hold off;
legend('angle','gamma');
%
%figure(2);
%plot(v,(pi-normal_angle)*180/pi,'b',v,(pi+normal_angle)*180/pi,'r');
%legend('pi-normal angle','pi+normal angle');
%
%figure(3);
%plot(v,d);
%legend('d');

