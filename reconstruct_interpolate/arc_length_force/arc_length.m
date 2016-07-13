clear all;
close all;
clc;

arc = [];
a = 1e-5;
while (a<0.2)
  arc = [arc;a];
  a = a/0.9;
end

f1 = 1./arc;

L = 0.3;
f2 = L-arc;

%b = 1e5;
b = 1e3;
f3 = b*arc;

ft = f1+f3;


%plot(arc,f1);
loglog(arc,ft,'*');
xlabel('arc length (um)');
ylabel('force');
hold on;
loglog(arc,f1,'ro');
%loglog(arc,f3,'ro');
%loglog(arc,ft,'ko');
hold off;
