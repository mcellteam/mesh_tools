clear all;
close all;
clc;

%[z,p,k] = buttap(n)
%[num,den] = zp2tf(z,p,k)
%bode(num,den)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% single copy of data set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mofo
%sampling_freq = 600;
%my_order = 3;
%cutoff_freq = 50 / (sampling_freq/2);
%[num,den] = butter(my_order,cutoff_freq,'low');
%%freqz(num,den,128,sampling_freq);
%smoothed = filtfilt(num,den,a);
%
%x = a(3);
%y = a(2);
%z = a(1);
%a = [x;y; z; a];
%mytry = [x y z];
%
%for i=4:length(a)
%  input = num(1)*a(i)+num(2)*a(i-1)+num(3)*a(i-2)+num(4)*a(i-3);
%  output = den(2)*mytry(i-1)+den(3)*mytry(i-2)+den(4)*mytry(i-3);
%  mytry = [mytry input-output];
%end
%
%b = fliplr(mytry(4:end));
%x = b(3);
%y = b(2);
%z = b(1);
%b = [z y x b];
%mytry2 = [z y x];
%for i=4:length(b)
%  input = num(1)*b(i)+num(2)*b(i-1)+num(3)*b(i-2)+num(4)*b(i-3);
%  output = den(2)*mytry2(i-1)+den(3)*mytry2(i-2)+den(4)*mytry2(i-3);
%  mytry2 = [mytry2 input-output];
%end
%
%final = fliplr(mytry2(4:end));
%
%figure(2);
%plot(a(4:end),'k+-')
%hold on;
%plot(smoothed,'ro-')
%plot(final,'bx')
%hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% multiple copies of data set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mofo3
sampling_freq = 600;
my_order = 3;
cutoff_freq = 50 / (sampling_freq/2);
[num,den] = butter(my_order,cutoff_freq,'low');
%freqz(num,den,128,sampling_freq);
d = [a' a' a']';
smoothed = filtfilt(num,den,d);

%x = d(3);
%y = d(2);
%z = d(1);
x = 0;
y = 0;
z = 0;
d = [x;y; z; d];
mytry = [x y z];

for i=4:length(d)
  input = num(1)*d(i)+num(2)*d(i-1)+num(3)*d(i-2)+num(4)*d(i-3);
  output = den(2)*mytry(i-1)+den(3)*mytry(i-2)+den(4)*mytry(i-3);
  mytry = [mytry input-output];
end

b = fliplr(mytry(4:end));
%x = b(3);
%y = b(2);
%z = b(1);
x = 0;
y = 0;
z = 0;
b = [z y x b];
mytry2 = [z y x];
for i=4:length(b)
  input = num(1)*b(i)+num(2)*b(i-1)+num(3)*b(i-2)+num(4)*b(i-3);
  output = den(2)*mytry2(i-1)+den(3)*mytry2(i-2)+den(4)*mytry2(i-3);
  mytry2 = [mytry2 input-output];
end

final = fliplr(mytry2(4:end));

j = 4+length(a);
k = j+length(a)-1;

newsmooth = circshift(smoothed(j:k),3);
newfinal = circshift(final(j:k)',3);

figure(1);
plot(d(j:k),'k+-')
hold on;
%plot(smoothed(j:k),'ro-')
plot(newsmooth,'ro-')
%plot(final(j:k),'bx')
%plot(newfinal,'g*')


%hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% gaussian convolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mofo2

%%%%%%%%%%%%%%% 3 elements

%c = [0.1587 0.6827 0.1587];
%
%o = [];
%for i=1:length(a)
%  if (i==1)
%    val = c(1)*a(end)+c(2)*a(i)+c(3)*a(i+1);
%  elseif (i==length(a))
%    val = c(1)*a(i-1)+c(2)*a(i)+c(3)*a(1);
%  else
%    val = c(1)*a(i-1)+c(2)*a(i)+c(3)*a(i+1);
%  end
%  o = [o val];
%end

%%%%%%%%%%%%%%% 5 elements

%c = [0.0668 0.2417 0.3829 0.2417 0.0668];
%b = [a(end-1);a(end);a;a(1);a(2)];
%
%o = [];
%for i=3:length(a)
%  val = c(1)*b(i-2)+c(2)*b(i-1)+c(3)*b(i)+c(4)*b(i+1)+c(5)*b(i+2);
%  o = [o val];
%end

%figure(3);
%plot(a,'k+-')
%hold on;
%plot(o,'ro-')
%hold off;

%%%%%%%%%%%%%%% 7 elements

%c = [0.1056 0.1210 0.1747 0.1974 0.1747 0.1210 0.1056];
%b = [a(end-2);a(end-1);a(end);a;a(1);a(2);a(3)];
%
%o = [];
%for i=4:length(a)
%  val = c(1)*b(i-3)+c(2)*b(i-2)+c(3)*b(i-1)+c(4)*b(i)+c(5)*b(i+1)+c(6)*b(i+2)+c(7)*b(i+3);
%  o = [o val];
%end
%
%figure(3);
%plot(a,'k+-')
%hold on;
%plot(o,'ro-')
%hold off;

%%%%%%%%%%%%%%% n elements

%n = 7;
%half_width = 1/(2^((n-1)/2-1));
%%[p0] = normspec([-half_width half_width]);
%p0 = normcdf(half_width,0,1)-normcdf(-half_width,0,1);
%
%m = (n-1)/2-1;
%p = [p0];
%for i=1:m,
%  g = 2*i-1;
%  h = g+2;
%  %[p1] = normspec([half_width*g half_width*h]);
%  p1 = normcdf(half_width*h,0,1)-normcdf(half_width*g,0,1);
%  p = [p1 p p1];
%end
%q = (1-sum(p))/2;
%p = [q p q];
%
%o = [];
%d = (n-1)/2;
%for i=1:length(a)
%  j = i-d-1;
%  val = 0;
%  for k=1:n,
%    index = j+k;
%    if (index<1)
%      index = index+length(a);
%    elseif (index>length(a))
%      index = index-length(a);
%    end
%    val = val+p(k)*a(index);
%  end
%  o = [o val];
%end

%close all;

%figure(2);
%plot(a,'k+-')
%hold on;
%plot(o,'g*-')
%hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% RADIUS OF CURVATURE STILL TOO CURVY %%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mofo2
sampling_freq = 600;
my_order = 3;
cutoff_freq = 75 / (sampling_freq/2);
[num,den] = butter(my_order,cutoff_freq,'low');
%freqz(num,den,128,sampling_freq);
d = [a' a' a']';
smoothed = filtfilt(num,den,d);

%x = d(3);
%y = d(2);
%z = d(1);
x = 0;
y = 0;
z = 0;
d = [x;y; z; d];
mytry = [x y z];

for i=4:length(d)
  input = num(1)*d(i)+num(2)*d(i-1)+num(3)*d(i-2)+num(4)*d(i-3);
  output = den(2)*mytry(i-1)+den(3)*mytry(i-2)+den(4)*mytry(i-3);
  mytry = [mytry input-output];
end

b = fliplr(mytry(4:end));
%x = b(3);
%y = b(2);
%z = b(1);
x = 0;
y = 0;
z = 0;
b = [z y x b];
mytry2 = [z y x];
for i=4:length(b)
  input = num(1)*b(i)+num(2)*b(i-1)+num(3)*b(i-2)+num(4)*b(i-3);
  output = den(2)*mytry2(i-1)+den(3)*mytry2(i-2)+den(4)*mytry2(i-3);
  mytry2 = [mytry2 input-output];
end

final = fliplr(mytry2(4:end));

j = 4+length(a);
k = j+length(a)-1;

newsmooth = circshift(smoothed(j:k),3);
newfinal = circshift(final(j:k)',3);

%figure(1);
%plot(d(j:k),'k+-')
%hold on;
%plot(smoothed(j:k),'ro-')
plot(newsmooth,'g*')
%plot(final(j:k),'bx')
%plot(newfinal,'g*')


hold off;





