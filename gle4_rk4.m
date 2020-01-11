% Stephen Wedekind
% 09 January 2019
% Generalized Lorenz System in 4D (GLE-4)
% Figure 1

%clear all; help gle4;
function [ BackPlot, FrontPlot, MasterPlot ] = gle4_rk4(x1, x2, x3, x4, t)

    % global beta rho gamma
    global timestep nstep transient_buffer

    % beta = 6;
    % rho = 8;
    % gamma = 0;

    % timestep = 0.01;
    % nstep = 230; %15000;

    %% 4th Order Runge-Kutta Integration:
    for i=1:nstep
     t(i+1) = t(i) + timestep;

     ka1 = x1dot( x1(i) , x2(i) , x3(i) , x4(i) , t(i) );
     ka2 = x1dot( ( (x1(i)) + (timestep/2)*(ka1) ) , ( (x2(i)) + (timestep/2)*(ka1) ) , ( (x3(i)) + (timestep/2)*(ka1) ) , ( (x4(i)) + (timestep/2)*(ka1) ) , ( (t(i)) + (timestep/2) ) );
     ka3 = x1dot( ( (x1(i)) + (timestep/2)*(ka2) ) , ( (x2(i)) + (timestep/2)*(ka2) ) , ( (x3(i)) + (timestep/2)*(ka2) ) , ( (x4(i)) + (timestep/2)*(ka2) ) , ( (t(i)) + (timestep/2) ) );
     ka4 = x1dot( ( (x1(i)) + (timestep)*(ka3) ) , ( (x2(i)) + (timestep)*(ka3) ), ( (x3(i)) + (timestep)*(ka3) ) , ( (x4(i)) + (timestep)*(ka3) ) , ( (t(i)) + (timestep) ) );

     kb1 = x2dot( x1(i) , x2(i) , x3(i) , x4(i) , t(i) );
     kb2 = x2dot( ( (x1(i)) + (timestep/2)*(kb1) ) , ( (x2(i)) + (timestep/2)*(kb1) ) , ( (x3(i)) + (timestep/2)*(kb1) ) , ( (x4(i)) + (timestep/2)*(kb1) ) , ( (t(i)) + (timestep/2) ) );
     kb3 = x2dot( ( (x1(i)) + (timestep/2)*(kb2) ) , ( (x2(i)) + (timestep/2)*(kb2) ) , ( (x3(i)) + (timestep/2)*(kb2) ) , ( (x4(i)) + (timestep/2)*(kb2) ) , ( (t(i)) + (timestep/2) ) );
     kb4 = x2dot( ( (x1(i)) + (timestep)*(kb3) ) , ( (x2(i)) + (timestep)*(kb3) ), ( (x3(i)) + (timestep)*(kb3) ) , ( (x4(i)) + (timestep)*(kb3) ) , ( (t(i)) + (timestep) ) );

     kc1 = x3dot( x1(i) , x2(i) , x3(i) , x4(i) , t(i) );
     kc2 = x3dot( ( (x1(i)) + (timestep/2)*(kc1) ) , ( (x2(i)) + (timestep/2)*(kc1) ) , ( (x3(i)) + (timestep/2)*(kc1) ) , ( (x4(i)) + (timestep/2)*(kc1) ) , ( (t(i)) + (timestep/2) ) );
     kc3 = x3dot( ( (x1(i)) + (timestep/2)*(kc2) ) , ( (x2(i)) + (timestep/2)*(kc2) ) , ( (x3(i)) + (timestep/2)*(kc2) ) , ( (x4(i)) + (timestep/2)*(kc2) ) , ( (t(i)) + (timestep/2) ) );
     kc4 = x3dot( ( (x1(i)) + (timestep)*(kc3) ) , ( (x2(i)) + (timestep)*(kc3) ), ( (x3(i)) + (timestep)*(kc3) ) , ( (x4(i)) + (timestep)*(kc3) ) , ( (t(i)) + (timestep) ) );

     kd1 = x4dot( x1(i) , x2(i) , x3(i) , x4(i) , t(i) );
     kd2 = x4dot( ( (x1(i)) + (timestep/2)*(kd1) ) , ( (x2(i)) + (timestep/2)*(kd1) ) , ( (x3(i)) + (timestep/2)*(kd1) ) , ( (x4(i)) + (timestep/2)*(kd1) ) , ( (t(i)) + (timestep/2) ) );
     kd3 = x4dot( ( (x1(i)) + (timestep/2)*(kd2) ) , ( (x2(i)) + (timestep/2)*(kd2) ) , ( (x3(i)) + (timestep/2)*(kd2) ) , ( (x4(i)) + (timestep/2)*(kd2) ) , ( (t(i)) + (timestep/2) ) );
     kd4 = x4dot( ( (x1(i)) + (timestep)*(kd3) ) , ( (x2(i)) + (timestep)*(kd3) ), ( (x3(i)) + (timestep)*(kd3) ) , ( (x4(i)) + (timestep)*(kd3) ) , ( (t(i)) + (timestep) ) ); 

     x1(i+1) = x1(i) + (timestep/6)*(ka1 + 2*ka2 + 2*ka3 + ka4);
     x2(i+1) = x2(i) + (timestep/6)*(kb1 + 2*kb2 + 2*kb3 + kb4);
     x3(i+1) = x3(i) + (timestep/6)*(kc1 + 2*kc2 + 2*kc3 + kc4);
     x4(i+1) = x4(i) + (timestep/6)*(kd1 + 2*kd2 + 2*kd3 + kd4);

     if ( i >= transient_buffer )
          x1plot(i) = x1(i+1);
          x2plot(i) = x2(i+1);
          x3plot(i) = x3(i+1);
          x4plot(i) = x4(i+1);
          if x1dot(x1(i+1), x2(i+1), x3(i+1), x4(i+1), t) < 0
              x1BackPlot(i) = x1(i+1);
              x2BackPlot(i) = x2(i+1);
              x3BackPlot(i) = x3(i+1);
              x4BackPlot(i) = x4(i+1);
          elseif x1dot(x1(i+1), x2(i+1), x3(i+1), x4(i+1), t) > 0
              x1FrontPlot(i) = x1(i+1);
              x2FrontPlot(i) = x2(i+1);
              x3FrontPlot(i) = x3(i+1);
              x4FrontPlot(i) = x4(i+1);
%           else
%               x1EdgePlot(i) = x1(i+1);
%               x2EdgePlot(i) = x2(i+1);
%               x3EdgePlot(i) = x3(i+1);
%               x4EdgePlot(i) = x4(i+1);
          end

     end

    end
    
    BackPlot = cell(1, 4);
    FrontPlot = cell(1, 4);
%     EdgePlot = cell(1, 4);
    MasterPlot = cell(1,4);
    
    BackPlot{1, 1} = x1BackPlot;
    BackPlot{1, 2} = x2BackPlot;
    BackPlot{1, 3} = x3BackPlot;
    BackPlot{1, 4} = x4BackPlot;
    
    FrontPlot{1, 1} = x1FrontPlot;
    FrontPlot{1, 2} = x2FrontPlot;
    FrontPlot{1, 3} = x3FrontPlot;
    FrontPlot{1, 4} = x4FrontPlot;
    
    MasterPlot{1, 1} = x1plot;
    MasterPlot{1, 2} = x2plot;
    MasterPlot{1, 3} = x3plot;
    MasterPlot{1, 4} = x4plot;
    
%     EdgePlot{1} = x1EdgePlot;
%     EdgePlot{2} = x2EdgePlot;
%     EdgePlot{3} = x3EdgePlot;
%     EdgePlot{4} = x4EdgePlot;

end
