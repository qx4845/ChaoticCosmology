% Stephen Wedekind
% 09 January 2019
% Generalized Lorenz System in 4D (GLE-4)
% Figure 1

%clear all; help gle4_fig1;

global beta
global rho
global gamma

beta = 6;%1;
rho = 8;%4.3;
gamma = 0;%1;

global timestep
global nstep

timestep = 0.01;
nstep = 300; %15000;

% x1 = NaN(nstep, 1);
% x2 = NaN(nstep, 1);
% x3 = NaN(nstep, 1);
% x4 = NaN(nstep, 1);
% t = NaN(nstep, 1);
% 
% x1(1) = 0.0;
% x2(1) = 1.0;
% x3(1) = 1.0;
% x4(1) = 1.0;
% t(1) = 0.0;

global transient_buffer nTries
transient_buffer = 0;%10000;

nTries = 10;

BackPlots = cell(nTries, 1);
FrontPlots = cell(nTries, 1);
% EdgePlots = cell(nTries, 1);
MasterPlots = cell(nTries, 1);

x1s = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
x2s = [0.1, 0.2, 0.3, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.0];
x3s = [0.6, 0.5, 0.4, 0.5, 0.2, 0.1, 0.2, 0.1, 0.0, 0.0, 0.5];
x4s = [0.3, 0.3, 0.3, 0.0, 0.3, 0.3, 0.1, 0.1, 0.1, 0.0, 0.5];

backCounter0 = 0;
frontCounter0 = 0;
masterCounter0 = 0;

for nTry=1:nTries
    
    x1 = NaN(nstep, 1);
    x2 = NaN(nstep, 1);
    x3 = NaN(nstep, 1);
    x4 = NaN(nstep, 1);
    t = NaN(nstep, 1);

    x1(1) = 0.0;
    x2(1) = x2s(nTry);
    x3(1) = x3s(nTry);
    x4(1) = x4s(nTry);
    t(1) = 0.0;
    
    [ BackPlots{nTry, 1}, FrontPlots{nTry, 1}, MasterPlots{nTry, 1} ] = gle4_rk4(x1, x2, x3, x4, t);
    
    backCounter0 = backCounter0 + length(BackPlots{nTry,1}{4});
    frontCounter0 = frontCounter0 + length(FrontPlots{nTry,1}{4});
    masterCounter0 = masterCounter0 + length(MasterPlots{nTry,1}{4});
end

MasterBackPlot = zeros(4, backCounter0);
MasterFrontPlot = zeros(4, frontCounter0);
MasterCombinedPlot = zeros(4, masterCounter0);

backCounter = 1;
frontCounter = 1;
masterCounter = 1;

for nTry=1:nTries
    
    MasterBackPlot(1, backCounter : backCounter + length(BackPlots{nTry,1}{1}) - 1) = BackPlots{nTry,1}{1};
    MasterBackPlot(2, backCounter : backCounter + length(BackPlots{nTry,1}{2}) - 1) = BackPlots{nTry,1}{2};
    MasterBackPlot(3, backCounter : backCounter + length(BackPlots{nTry,1}{3}) - 1) = BackPlots{nTry,1}{3};
    MasterBackPlot(4, backCounter : backCounter + length(BackPlots{nTry,1}{4}) - 1) = BackPlots{nTry,1}{4};
    
    MasterFrontPlot(1, frontCounter : frontCounter + length(FrontPlots{nTry,1}{1}) - 1) = FrontPlots{nTry,1}{1};
    MasterFrontPlot(2, frontCounter : frontCounter + length(FrontPlots{nTry,1}{2}) - 1) = FrontPlots{nTry,1}{2};
    MasterFrontPlot(3, frontCounter : frontCounter + length(FrontPlots{nTry,1}{3}) - 1) = FrontPlots{nTry,1}{3};
    MasterFrontPlot(4, frontCounter : frontCounter + length(FrontPlots{nTry,1}{4}) - 1) = FrontPlots{nTry,1}{4};
    
    MasterCombinedPlot(1, masterCounter : masterCounter + length(MasterPlots{nTry,1}{1}) - 1) = MasterPlots{nTry,1}{1};
    MasterCombinedPlot(2, masterCounter : masterCounter + length(MasterPlots{nTry,1}{2}) - 1) = MasterPlots{nTry,1}{2};
    MasterCombinedPlot(3, masterCounter : masterCounter + length(MasterPlots{nTry,1}{3}) - 1) = MasterPlots{nTry,1}{3};
    MasterCombinedPlot(4, masterCounter : masterCounter + length(MasterPlots{nTry,1}{4}) - 1) = MasterPlots{nTry,1}{4};
    
    backCounter = backCounter + length(BackPlots{nTry,1}{4});
    frontCounter = frontCounter + length(FrontPlots{nTry,1}{4});
    masterCounter = masterCounter + length(MasterPlots{nTry,1}{4});
    
end

figure (1); clf;
% plot(x2plot,x3plot,'.');
plot(MasterFrontPlot(2,:),MasterFrontPlot(4,:),'r.');
hold on;
plot(MasterBackPlot(2,:),MasterBackPlot(4,:),'b.');
xlabel('x2 - position');
ylabel('x4 - position');
title('x2x4 - projection (x4 vs. x2)');
%text(0,-1,"x0 = 2.5 , F = 7.5");

figure (2); clf;
% plot(x2plot,x3plot,'.');
plot(MasterFrontPlot(1,:),MasterFrontPlot(2,:),'r.');
hold on;
plot(MasterBackPlot(1,:),MasterBackPlot(2,:),'b.');
xlabel('x2 - position');
ylabel('x1 - position');
title('x1x2 - projection (x2 vs. x1)');

% 
pause(1);
refresh(1);
figure(3);  clf;
% plot3(x1plot,x2plot,x3plot);  axis([-3,3,-3,3,-3,3]);
% plot3(MasterCombinedPlot(1,:),MasterCombinedPlot(2,:),MasterCombinedPlot(3,:));  axis([-3,3,-3,3,-3,3]);
plot3(MasterFrontPlot(1,:),MasterFrontPlot(2,:),MasterFrontPlot(3,:), 'r');  %axis([-3,3,-3,3,-3,3]);
hold on;
plot3(MasterBackPlot(1,:),MasterBackPlot(2,:),MasterBackPlot(3,:), 'b');  %axis([-3,3,-3,3,-3,3]);
xlabel('x');  ylabel('y');  zlabel('z');
title('GLE-4 (3D Projection)');
refresh(1);