function [] = plot_optimal_control_and_variables_over_time(inputs)
%plot_optimal_control_and_variables_over_time Creates a plot of the output
%from the optimal control calculation -
% Plots optimal control and other variables
%over time. 
% copied from dara_mm_initial_sims.m
% Originally adapted from the code for fig 7
% this version designed to be multi-purpose
% takes dictionary of inputs, plot depends on what is provided and flag
% values; can add options to extend with backwards compat

% set any defaults, then dump the inputs map into local namespace,
% overriding defaults if a value is provided
plotPNsum = true;
includeLegend = true;
xLab = "Time";
yLab = "State";
for key = inputs.keys()
    eval(append(key{1}," = inputs('",key{1},"');"));
end

colours = [ 
    30/255  136/255  229/255 % blue; cost (and A)
    216/255  27/255  96/255 % red; P
    255/255  193/255  7/255 % yellow; N  
];

%colours = [ 
%    0/255  114/255  189/255
%    222/255  125/255  0/255 
%    220/255  220/255  10/255
%]; %Define colours for plot
if plotPNsum
    colours = [colours ; 237/255  177/255  32/255 ];
end

% inputs = containers.Map({'a','b'},{5,3:6});
if isKey(inputs,"saveName")
    figure('Name',saveName );
else
    figure;
end
set(gca, 'ColorOrder', colours);
hold on
box on
lineA = plot(t_y,y(1,:),'LineWidth',2);
lineP = plot(t_y,y(2,:),'LineWidth',2);
lineN = plot(t_y,y(3,:),'LineWidth',2);
lineArr = [lineA,lineP,lineN];
lineLabs = {'A','P','N'};
if plotPNsum
    linePNsum = plot(t_y,y(2,:)+y(3,:),'LineWidth',2);
    lineArr(end+1) = linePNsum;
    lineLabs{end+1} = 'P+N';
end
lineU = plot(t_y,U,'k','LineWidth',2);
lineArr(end+1) = lineU;
lineLabs{end+1} = 'u';
% hL = legend([line1,line2,line3,line4,line5,line6],{'A','P','N','P+N','u','psi'},'Location','northeast');
if includeLegend
    legend(lineArr,lineLabs,'Location','northeast');
end
ylabel(yLab,'fontsize',18);
xlabel(xLab,'fontsize',18);
axis([0,t_y(end),0,1]) % axis([0,Tfinal,0,1])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)

end