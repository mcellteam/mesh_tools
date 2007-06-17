clear all;
close all;
clc;

% number of samples
n = 100;
v = [1:n];

% angle of edge, radians
theta = linspace(0,2*pi,n)';

v1 = [1 1 1];
v2 = [1 1 10];
o1 = [3 1 3];
o2 = [2*cos(theta)+1 2*sin(theta)+1 7*ones(n,1)];

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
phi = acos(dot(repmat(s1,[n,1]),s2,2)./sqrt(dot(repmat(s1,[n,1]),repmat(s1,[n,1]),2))./sqrt(dot(s2,s2,2)));

gamma = [gamma_first(1:50);-gamma_first(51:100)];

% plots
figure(1);
plot(v,theta*180/pi,'k',v,phi*180/pi,'b',v,gamma*180/pi,'r',v,(gamma+theta)*180/pi,'c');
legend('theta','phi','gamma','gamma+theta');
