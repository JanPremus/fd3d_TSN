%Print time series of traction and slip rate in desired point on the fault
%Beggining at the left down corner 
nx = 20; %Horizontal position on the fault in km 
nz = 8;  %Vertical position on the fault in km 
f = fopen('inputfd3d.dat');
str = fgetl(f);
size(1:3) = sscanf(str,'%d',3);
str = fgetl(f);
dh = sscanf(str,'%f',1);
str = fgetl(f);
NT = sscanf(str,'%f',1);
str = fgetl(f);
dt = sscanf(str,'%f',1);
fclose(f);

PX=floor(nx*1000/dh);
PZ=floor(nz*1000/dh);

f(1) = fopen('result/shearstressX.res');
f(2) = fopen('result/shearstressZ.res');
f(3) = fopen('result/sliprateX.res');
f(4) = fopen('result/sliprateZ.res');

for j=1:4
    for i=1:NT
        D = fread(f(j),size(1)*size(3),'real*8');
        D = reshape(D(1:size(1)*size(3)),[size(1) size(3)]);
        series{j}(i)=D(PX,PZ);
    end
end

time = dt*(1:NT);




figure

sp1=subplot(4,1,1);
plot(time,series{1});
title('Horizontal traction (Pa)')
pbaspect([2 1 1])

sp2=subplot(4,1,2);
plot(time,series{3});
title('Horizontal slip rate (m/s)')
pbaspect([2 1 1])

sp3=subplot(4,1,3);
plot(time,series{2});
title('Vertical traction (Pa)')
pbaspect([2 1 1])

sp4=subplot(4,1,4);
plot(time,series{4});
title('Vertical slip rate (m/s)')
pbaspect([2 1 1])

print('TimeSeries','-dpng')

