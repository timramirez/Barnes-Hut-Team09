clear all;
close all;

nBod = 3;
data = [6.10357e+00  3.83627e+00  0.00000e+00  0.00000e+00  8.04789e+01
 6.00033e+00  3.99895e+00  0.00000e+00  0.00000e+00  2.77796e+01
 5.99765e+00  4.00397e+00  0.00000e+00  0.00000e+00  3.49997e+01];
pos = data(:,1:2);
mass = data(:,end);
F = [0 0; 0 0; 0 0];
for iBod = 1:nBod
    if iBod < nBod
        for jBod = (iBod+1):nBod
            R = sqrt((pos(iBod, 1)-pos(jBod, 1))^2 + (pos(iBod, 2)-pos(jBod, 2))^2);
            uij = (pos(jBod,:)-pos(iBod,:))/R;
            Fij = (mass(iBod)*mass(jBod)/(R^2)).*uij;
            F(iBod,:) = F(iBod,:)+Fij;
            F(jBod,:) = F(jBod,:)-Fij;
        end
    end
end