function [timespan, species, kinaseSignals] = run_model()
% Function to run AP1 regulatory model. All the dependent or nested 
% functions are in this file
%
% To run the simulation, type 
% >> [t, S, K] = run_model; 
%
% Output: timespan      - time vector
%         species       - species object that has
%                         .names - names of the species
%                         .IC - initial condition
%                         .steadyStateC - steady state condition
%                         .activityLevels - kinase stimulation
%         kinaseSignals - kinase signals
%
%


format long eng

options = odeset( 'RelTol', 1e-9, 'AbsTol', 1e-9 );
timespan = 0:0.1:60; % regulatory time span of one hour

%% % Initialization of kinase feature values and initial concentration of species

% kinase features and their values

kinase.featureNames={'1, a1F'
'2, a2F'
'3, a3F'
'4, t1F'
'5, t2F'
'6, t3F'
%
'7, a1E'
'8, a2E'
'9, a3E'                 
'10, t1E'
'11, t2E'
'12, t3E'
%
'13, a1J'
'14, a2J'
'15, a3J'
'16, t1J'
'17, t2J'
'18, t3J'
};

kinase.featuresValues = [
2
13
15
3
130
20
0
8
20
1.5
17.5
5000
5
25
15
7
22
20
];

% species names

species.names={'1, cFos_nucleus'
'2, pcFOS_nucleus'
'3, ppcFOS_nucleus'
%
'4, ELK-1_nucleus' %y0[4]  = 56.3120
'5, pELK-1_nucleus'
%
'6, cJUN_nucleus'
'7, pcJUN_nucleus'
'8, ppcJUN_nucleus'
%
'9, ATF-2_nucleus'  %y0[9]  = 28.1560                 
'10, pATF-2_nucleus'
'11, ppATF-2_nucleus'
%
'12, (ppcFos:ppcJun)_nucleus'
'13, (ppcJun:ppcJun)_nucleus'
'14, (ppcJun:ppATF-2)_nucleus'
%
'15, cFos(promoter)_nucleus' %y0[15] = 0.2350
'16, (cFos(promoter):pELK-1)_nucleus'
'17, (cFos(pre)-mRNA)_nucleus'
'18, (cFos-mRNA)_cytosol'
'19, cFOS_cytosol'
%
'20, cJun(promoter)_nucleus' %y0[20] = 0.2350
'21, (cJun(promoter):cJun:ATF-2)_nucleus'
'22, (cJun(pre)-mRNA)_nucleus'
'23, (cJun-mRNA)_cytosol'
'24, cJUN_cytosol'
%
'25, DownstreamGenes(promoter)_nucleus' %y0[25] = 0.2350
'26, (DownstreamGenes(promoter):cJun:cJun)_nucleus'
'27, (DownstreamGenes(promoter):cFos:cJun)_nucleus'
'28, (DownstreamGenes(pre)-mRNA)_nucleus'
'29, (DownstreamGenes-mRNA)_cytosol'
%
'30, Total_AP-I = (cFos:cJun)_nucleus + (cJun:cJun)_nucleus'
};

% initial condition
species.IC     = zeros(1,30);
species.IC(4)  = 56.3120;
species.IC(9)  = 28.1560;
species.IC(15) = 0.2350;
species.IC(20) = 0.2350;
species.IC(25) = 0.2350;

%% % Generate signaling kinases

kinaseSignals = genKinaseFun(kinase.featuresValues);
kinaseSignals = [kinaseSignals(:,1), 100*kinaseSignals(:,2), kinaseSignals(:,1), 100*kinaseSignals(:,3), kinaseSignals(:,1), 100*kinaseSignals(:,4)];
kinaseSignals = downsample(kinaseSignals,10);

%% % Run simulation
[ ~, species.steadyStateC ] = ode15s( @AP1_reg_model, [0 10^4], species.IC, options);
[ timespan, species.activityLevels ] = ode15s( @AP1_reg_model,[0 60], species.steadyStateC(end,:), options, kinaseSignals );

%% Generate plots
Linept = 3;
subplot(2,4,1)
l1 = plot(kinaseSignals(:,3), kinaseSignals(:,4), 'LineWidth', Linept);
title('ERK');
ylabel('Activity levels');

subplot(2,4,2)
plot(kinaseSignals(:,1), kinaseSignals(:,2), 'LineWidth', Linept);
title('FRK');

subplot(2,4,3)
plot(kinaseSignals(:,5), kinaseSignals(:,6), 'LineWidth', Linept);
title('JNK');

subplot(2,4,4)
l2 = plot(timespan, species.activityLevels(:,30), 'r', 'LineWidth', Linept);
title('Total AP-1');

subplot(2,4,5)
plot(timespan, species.activityLevels(:,18), 'r', 'LineWidth', Linept);
title('c-Fos mRNA');
xlabel('Time (min)');
ylabel('Activity levels');

subplot(2,4,6)
plot(timespan, species.activityLevels(:,23), 'r', 'LineWidth', Linept);
title('c-Jun mRNA');
xlabel('Time (min)');

subplot(2,4,7)
plot(timespan, species.activityLevels(:,12), 'r', 'LineWidth', Linept);
title('Heterodimer');
xlabel('Time (min)');

subplot(2,4,8)
plot(timespan, species.activityLevels(:,13), 'r', 'LineWidth', Linept);
title('Homodimer');
xlabel('Time (min)');

hL = legend([l1,l2],{'Input Signals','Downstream response'},'Orientation','horizontal');
newPosition = [0.4 0.4 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
end



function kinaseSignals = genKinaseFun( kinaseFeatures )
% Phenomological model to generate kinase signals
%
% Output: kinaseSignals - kinase activity levels as time series
%

kinaseSignals = zeros(601,4);
kinaseSignals(:,1) = 0:0.1:60;

  for i=0:2

    j1 = ceil( 10 * kinaseFeatures(1 + i*6));
    j2 = j1 +  ceil( 10 * kinaseFeatures(2 + i*6));
    j3 = j2 +  ceil( 10 * kinaseFeatures(3 + i*6));

    kinaseSignals((j1+1):j2,2+i) = 1 - exp( -( kinaseSignals((j1+1):j2,1) - kinaseFeatures(1+i*6) ) / kinaseFeatures(4 + i*6) );
    kinaseSignals((j2+1):j3,2+i) = kinaseSignals(j2,2+i) * exp( -( kinaseSignals((j2+1):j3,1) - kinaseFeatures(1+i*6) - kinaseFeatures(2 + i*6) ) / kinaseFeatures( 5 + i*6 ) );
    kinaseSignals((j3+1):601,2+i) = kinaseSignals(j3,2+i) * exp( -( kinaseSignals((j3+1):601,1) - kinaseFeatures(1+i*6) - kinaseFeatures(2 + i*6) - kinaseFeatures(3 + i*6) ) / kinaseFeatures(6 + i*6) );

  end

end
 

