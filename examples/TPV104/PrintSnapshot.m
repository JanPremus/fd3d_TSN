%Print snapshots of traction and slip rate in desired time
T  = 3.123; %Time of the snapshot

f = fopen('inputfd3d.dat');
str = fgetl(f);
size(1:3) = sscanf(str,'%d',3);
str = fgetl(f);
str = fgetl(f);
str = fgetl(f);
dt = sscanf(str,'%f',1);
fclose(f);

f(1) = fopen('result\shearstressX.res');
f(2) = fopen('result\shearstressZ.res');
f(3) = fopen('result\sliprateX.res');
f(4) = fopen('result\sliprateZ.res');

TL=floor(T/dt);
for j=1:4
    for i=1:TL
        D{j} = fread(f(j),size(1)*size(3),'real*4');
        D{j} = reshape(D{j}(1:size(1)*size(3)),[size(1) size(3)]);

    end
end

figure

sp1=subplot(4,1,1);
imagesc(rot90(D{1},1));
title('Horizontal traction (Pa)')
colorbar;
pbaspect([size(1) size(3) 1])

sp2=subplot(4,1,2);
imagesc(rot90(D{3},1));
title('Horizontal slip rate (m/s)')
colorbar;
pbaspect([size(1) size(3) 1])

sp3=subplot(4,1,3);
imagesc(rot90(D{2},1));
title('Vertical traction (Pa)')
colorbar;
pbaspect([size(1) size(3) 1])

sp4=subplot(4,1,4);
imagesc(rot90(D{4},1));
title('Vertical slip rate (m/s)')
colorbar;
pbaspect([size(1) size(3) 1])

print('Snapshot','-dpng')
