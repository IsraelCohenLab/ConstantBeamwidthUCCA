% this script is used to extract all the figures presented in the paper
% "Constant-Beamwidth Beamforming with Concentric Ring Arrays", published
% in Sensor. Each section in the script extracts different graph.
close all; clear all; clc
%% Load pre-saved data
LP_UCCA_data = load('saved_results\Low_Pass_UCCA');             % Original constant BW design       
BP_UCCA_data = load('saved_results\Modified_UCCA');             % Low frequencies modification
Imp_UCCA_data = load('saved_results\Improved_Modified_UCCA');   % DI improvement and low frequencies modification
ULA_data = load('saved_results\ULA');                           % Constant beamwidth ULA
close all;
%% Plot beampatterns
f = LP_UCCA_data.frequency_vec;     % saved frequency axis
phi = (-90:1:90)';
% Plot beampattern of final result Imp_UCCA
figure;
imagesc(phi,f,20*log10(Imp_UCCA_data.Da'));
axis 'xy'

ax = gca;
new_label = arrayfun( @(x) num2str(x/1000) , ax.YTick, 'UniformOutput', false );
ax.YTickLabel = new_label;

colorbar;
colormap jet
caxis([-40 0]) 

new_label = arrayfun( @(x) num2str(x/1000) , ax.YTick, 'UniformOutput', false );
ax.YTickLabel = new_label;
v = [-3,-3];
hold on;
contour(phi,f,20*log10(LP_UCCA_data.Da'),v,'LineColor','white','LineWidth',4,'LineStyle',':')
% Beampattern compare of low frequency modification
% ax.YLim = [0 3000];
figure;
subplot(1,2,1)
imagesc(phi,f,20*log10(LP_UCCA_data.Da'));
axis 'xy'
ax = gca;
new_label = arrayfun( @(x) num2str(x/1000) , ax.YTick, 'UniformOutput', false );
ax.YTickLabel = new_label;
colormap jet
caxis([-40 0]) 
hold on;
contour(phi,f,20*log10(LP_UCCA_data.Da'),v,'LineColor','white','LineWidth',4,'LineStyle',':')
 ax.YLim = [0 3000];
 
subplot(1,2,2)
imagesc(phi,f,20*log10(BP_UCCA_data.Da'));
axis 'xy'
ax = gca;
new_label = arrayfun( @(x) num2str(x/1000) , ax.YTick, 'UniformOutput', false );
ax.YTickLabel = new_label;
colormap jet
caxis([-40 0]) 
hold on;
contour(phi,f,20*log10(BP_UCCA_data.Da'),v,'LineColor','white','LineWidth',4,'LineStyle',':')
 ax.YLim = [0 3000];
%% Define preformance graph parameters
linewd = 1;
MarkerSize=8;
hcfontsize = 12;
position1=[0.1300    0.7093    0.4    0.2157]';
position2=[0.1300    0.4596    0.4    0.2157]';
color = [0    0.5    0.0];
%% Extract DF and WNG graphs
fig = figure;
fig.InnerPosition = [680 228 843 750];
fig.OuterPosition = [680 228 843 750];
%%%%%%%% Plot DI of results %%%%%%%%
subplot(2,1,1)
ax1 = plot(f/1e3,10*log10(LP_UCCA_data.total_DF),'-','linewidth',linewd, 'Color',color);
hold on
ax2 = plot(f/1e3,10*log10(BP_UCCA_data.total_DF),'--b','linewidth',linewd);
ax3 = plot(f/1e3,10*log10(Imp_UCCA_data.total_DF),':r','linewidth',2, 'MarkerSize',MarkerSize);
ax4 = plot(f/1e3,10*log10(ULA_data.total_DF),'-.m','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
% set figure parameters
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(fig.Children, 'XScale', 'log')
set(gca,'XTickLabel',[])
box on; grid on;
% add markers on log-scaled axis for each of the graphs
x_mark = logspace(log10(ax1.XData(2)),log10(ax1.XData(end)),25);
y_mark = interp1(ax1.XData(2:end),ax1.YData(2:end),x_mark);
f_mark = figure(500);
figure(500); 
plot(x_mark, y_mark, 'o','MarkerSize',MarkerSize,'Color',color )
set(f_mark.Children, 'XScale', 'log')
copyobj(f_mark.Children.Children, fig.Children)

x_mark = logspace(log10(ax2.XData(2)),log10(ax2.XData(end)),25);
y_mark = interp1(ax2.XData(2:end),ax2.YData(2:end),x_mark);
f_mark = figure(500);
figure(500); 
plot(x_mark, y_mark, 'b*','MarkerSize',MarkerSize )
set(f_mark.Children, 'XScale', 'log')
copyobj(f_mark.Children.Children, fig.Children)

x_mark = logspace(log10(ax3.XData(2)),log10(ax3.XData(end)),25);
y_mark = interp1(ax3.XData(2:end),ax3.YData(2:end),x_mark);
f_mark = figure(500);
figure(500); 
plot(x_mark, y_mark, 'rs','MarkerSize',MarkerSize )
set(f_mark.Children, 'XScale', 'log')
copyobj(f_mark.Children.Children, fig.Children)

x_mark = logspace(log10(ax4.XData(2)),log10(ax4.XData(end)),25);
y_mark = interp1(ax4.XData(2:end),ax4.YData(2:end),x_mark);
f_mark = figure(500);
figure(500); 
plot(x_mark, y_mark, 'm^','MarkerSize',MarkerSize )
set(f_mark.Children, 'XScale', 'log')
copyobj(f_mark.Children.Children, fig.Children)

close(f_mark)
%%%%%%%% Plot WNG of results %%%%%%%%
sub_fig = subplot(2,1,2);
hold on
ax1 = plot(f/1e3,10*log10(LP_UCCA_data.total_DF_WNG),'-','linewidth',linewd,'Color',color);
ax2 = plot(f/1e3,10*log10(BP_UCCA_data.total_DF_WNG),'--b','linewidth',linewd);
ax3 = plot(f/1e3,10*log10(Imp_UCCA_data.total_DF_WNG),':r','linewidth',2);
ax4 = plot(f/1e3,10*log10(ULA_data.total_DF_WNG),'-.m','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off

set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca, 'YLim', [0 22]); 
set(gca, 'XScale', 'log')
fig.Children(1).XTick = [0 1 8];
set(gca,'XTickLabel',{'0','1','8'})
box on; grid on;

x_mark = logspace(log10(ax1.XData(2)),log10(ax1.XData(end)),25);
y_mark = interp1(ax1.XData(2:end),ax1.YData(2:end),x_mark);
f_mark = figure(500);
figure(500); 
plot(x_mark, y_mark, 'o','MarkerSize',MarkerSize,'Color',color )
set(f_mark.Children, 'XScale', 'log')
copyobj(f_mark.Children.Children, fig.Children(1))

x_mark = logspace(log10(ax2.XData(2)),log10(ax2.XData(end)),25);
y_mark = interp1(ax2.XData(2:end),ax2.YData(2:end),x_mark);
f_mark = figure(500);
figure(500); 
plot(x_mark, y_mark, 'b*','MarkerSize',MarkerSize )
set(f_mark.Children, 'XScale', 'log')
copyobj(f_mark.Children.Children, fig.Children(1))


x_mark = logspace(log10(ax3.XData(2)),log10(ax3.XData(end)),25);
y_mark = interp1(ax3.XData(2:end),ax3.YData(2:end),x_mark);
f_mark = figure(500);
figure(500); 
plot(x_mark, y_mark, 'rs','MarkerSize',MarkerSize )
set(f_mark.Children, 'XScale', 'log')
copyobj(f_mark.Children.Children, fig.Children(1))


x_mark = logspace(log10(ax4.XData(2)),log10(ax4.XData(end)),25);
y_mark = interp1(ax4.XData(2:end),ax4.YData(2:end),x_mark);
f_mark = figure(500);
figure(500); 
plot(x_mark, y_mark, 'm^','MarkerSize',MarkerSize )
set(f_mark.Children, 'XScale', 'log')
copyobj(f_mark.Children.Children, fig.Children(1))
close(f_mark)
%% Extract beamwidth comparison graph
phi = (-90:1:90)';
desired_bw_point = BP_UCCA_data.desired_bw_point;
total_data = {LP_UCCA_data, BP_UCCA_data, Imp_UCCA_data, ULA_data};
% define graph characteraristics
colors = {color, 'b' , 'r','m'};
markers = {'o','*','s','^'};
line_styles = {'-','--',':','-.'};
line_width = 1;

fig = figure; 
hold on
for i = 1:4
    curr_data = total_data{i};
    % unpack current plot charectaristics
    line_style = line_styles{i};
    marker = markers{i};
    line_color = colors{i};
    % adjust linewidth for visualization 
    if i ==3
        line_width = 2;
    else
        line_width = 1;
    end
    % calculate the effective beamwidth based on the minimal angle forwhich
    % the beampattern amplitude is closest to the desired
    difference_mat = abs(curr_data.Da - desired_bw_point);
    [min_points_values,min_points_id] = min(difference_mat,[],1);
    curr_ax = plot(f/1e3,2*abs(phi(min_points_id)),'LineStyle',line_style,'Color',line_color,'LineWidth',line_width);
    % add markers to the current graph
    x_mark = logspace(log10(curr_ax.XData(2)),log10(curr_ax.XData(end)),25);
    y_mark = interp1(curr_ax.XData(2:end),curr_ax.YData(2:end),x_mark);
    f_mark = figure(500);
    figure(500); 
    plot(x_mark, y_mark, 'Color', line_color, 'Marker', marker ,'MarkerSize',MarkerSize ,'LineStyle','none')
    set(f_mark.Children, 'XScale', 'log')
    copyobj(f_mark.Children.Children, fig.Children)
    
    close(f_mark)
end
hold off

set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca, 'XLim', [0 f(end-2)/1e3]); 
set(gca, 'XScale', 'log')
fig.Children(1).XTick = [0 1 f(end-2)/1e3];
set(gca,'XTickLabel',{'0','1','8'})
box on; grid on;
%% Sidelobe attenuation graph
DA_ULA = ULA_data.Da';
DA_CCA = Imp_UCCA_data.Da';
theta_vec = -90:90;

min_y = min(min(20*log10(DA_ULA)));
figure;
cc = jet(10);
BW = 15;
% plot beampattern of several frequencies
for id = 45:50
    subplot(1,2,1)
    hold on
    plot(theta_vec,20*log10(DA_CCA(id,:)),'b','LineWidth',0.5)
    plot([BW BW],[min_y 0],'--k')
    plot([-BW -BW],[min_y 0],'--k')
    plot(theta_vec,-3*ones(length(theta_vec),1),'--g','LineWidth',1)
    plot(theta_vec,-18*ones(length(theta_vec),1),'--m','LineWidth',2)

    xlabel('\theta [deg]')
    ylabel('Gain [dB]')
    axis tight
    ylim([min_y 0])
    
    subplot(1,2,2)
    hold on
    plot(theta_vec,20*log10(DA_ULA(id,:)),'r','LineWidth',0.5)
    plot([BW BW],[min_y 0],'--k')
    plot([-BW -BW],[min_y 0],'--k')
    plot(theta_vec,-3*ones(length(theta_vec),1),'--g','LineWidth',1)
    plot(theta_vec,-11*ones(length(theta_vec),1),'--m','LineWidth',2)


    drawnow
end
hold off
xlabel('\theta [deg]')
ylabel('Gain [dB]')
axis tight
ylim([min_y 0])
%% Compare Design compatability to various beampattern
linewd = 1.5;
MarkerSize=8;
hcfontsize = 12;

for i = 20:10:50
    load(['saved_results/' ,num2str(i),'.mat'])
    phi = (-90:1:90)';
    Vq = interp2(Da,5);
    figure;
    imagesc(phi,freq_vec,20*log10(Vq'));
    axis 'xy'
    ax = gca;
    new_label = arrayfun( @(x) num2str(x/1000) , ax.YTick, 'UniformOutput', false );
    ax.YTickLabel = new_label;
    colormap jet
    caxis([-40 0])
    v = [-3,-3];
    hold on;
    contour(phi,freq_vec,20*log10(Da'),v,'LineColor','white','LineWidth',4,'LineStyle',':')
    hold off;
    title(['BW = ',num2str(i)])
    set(gca, 'Color', [1, 1, 1]); 
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', hcfontsize);
    set(gca, 'LineWidth', linewd);
end

%% Plot DI and WNG of different beamwidths
BW_vec = 20:10:50;
Line_styles = {'-','-.','--',':'};
%Colors = []

fig = figure;
fig.InnerPosition = [680 228 843 750];
fig.OuterPosition = [680 228 843 750];
subplot(2,1,1)
hold on
for i = 1:4
    load(['saved_results\' ,num2str(BW_vec(i)),'.mat'])
    ax = plot(freq_vec/1e3 , total_DF,'linewidth',2,'LineStyle',Line_styles{i});
end
hold off

gca = fig.Children;
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', 1); 
 set(fig.Children, 'XScale', 'log')
box on; grid on;
set(fig.Children,'XTickLabel',[])
set(fig.Children,'XLim',[0,8])

sub_fig = subplot(2,1,2);
hold on
for i = 1:4
    load(['saved_results\' ,num2str(BW_vec(i)),'.mat'])
    ax = plot(freq_vec/1e3 , total_DF_WNG,'linewidth',2,'LineStyle',Line_styles{i});
end
hold off

gca = fig.Children(1);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', 1); 
set(gca, 'XScale', 'log')

fig.Children(1).XTick = [0 1 8];
set(fig.Children(1),'XTickLabel',{'0','1','8'})
box on; grid on;
set(fig.Children(1),'XLim',[0,8])