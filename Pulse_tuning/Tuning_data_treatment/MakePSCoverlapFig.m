function [] = MakePSCoverlapFig()
close all

load PulseData.mat

A = PulseD(:,1)
B = PulseD(:,2)

figure('color','w')

x = [0.2 3.84+0.5 3.84+0.5 1+0.5];
x = [-0.5 0.5 0.4 -0.4]+0.3;

cut = 23;
xx = [A(1:cut); 0.8; 0]
yy = [B(1:cut); 0; 0]


y = [0 0 12.72 12.72];
patchColor = [.8 .8 .8];
p1 = patch(x,y,patchColor)

patchColor2 = [.6 .6 .6];
p2 = patch(xx,yy,patchColor2)
%line([0.3 0.3], [0.0 12.72],'linestyle', ':', 'linewidth', 0.75, 'color', 'k')
hold on

p3 = plot(A, B, 'color', 'b', 'linewidth', 2)

xlabel('Time [ms]', 'fontname', 'times', 'fontsize', 18)
ylabel('Intensity [a.u.]', 'fontname', 'times','fontsize', 18)
set(gca, 'xlim', [-1 6])
legend([p1 p2 p3], 'Nominal PSC opening', 'Realized PSC opening', 'Moderator pulse') 
StandardFigureStuff(16)
