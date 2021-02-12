%min= load('max-min.dat');
%close;plot(min)

i = 1:length(min);
i1 = 1:length(min1);
i2 = 1:length(min2);
i3 = 1:length(min3);
close all; plot(i, min, i1, min1, i2, min2, i3, min3);
axis([0 10000 0.85 1.05]);
xlabel('iteration');
ylabel('min');
legend('Np=1', 'Np=2', 'Np=4', 'Np=8');
