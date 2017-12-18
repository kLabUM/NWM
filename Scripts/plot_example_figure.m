clear all
close all
cd('~/Google Drive/Docs Kevin/National Water Model/NWM/Scripts');
load('mar29.mat');
data = mar29;
clear mar29
%%
% bad_ind = find(data.ifis_id=='RAPIDCR01');
% bad_ind = find(data.ifis_id=='LIMECR01');
bad_ind = find(data.ifis_id=='WLNTCRK01');
% bad_ind = find(data.ifis_id=='BEAVER02');
good_ind = find(data.ifis_id=='CLRCRK01');

subplot(2,2,1)
scatter(data.t0(good_ind),data.gage_height(good_ind));
xlabel('Flow (cfs)')
ylabel('Water Level (ft)')
% title('Clear Creek');
set(gca,'FontSize',16)

subplot(2,2,2)
scatter(data.t0(bad_ind),data.gage_height(bad_ind));
xlabel('Flow (cfs)')
ylabel('Water Level (ft)')
% title('Walnut Creek');
set(gca,'FontSize',16)

subplot(2,2,3)
yyaxis left
plot(data.gage_height(good_ind),'LineWidth',2)
ylabel('Water Level (ft)')
yyaxis right
plot(data.t0(good_ind),'LineWidth',2)
ylabel('Flow (cfs)')
legend('Measured Height','Modeled Flow')
xlabel('Time (Hrs)')
set(gca,'FontSize',16)

subplot(2,2,4)
yyaxis left
plot(data.gage_height(bad_ind),'LineWidth',2)
ylabel('Water Level (ft)')
yyaxis right
plot(data.t0(bad_ind),'LineWidth',2)
ylabel('Flow (cfs)')
legend('Measured Height','Modeled Flow')
xlabel('Time (Hrs)')
set(gca,'FontSize',16)
%%
clear all, close all
load ss_data_TF_only.mat %change filename of data from first section as appropriate
data.ifis_id = string(data.ifis_id);
ifis_nwm = readtable('IFIS_NWM_data_culled.txt');
k = 1;
omit_ind = [];
for i = 1:length(data.ifis_id)
    ind = find(ifis_nwm.ifis_id==data.ifis_id(i));
    if length(ind)>0
        nwm_id(k) = ifis_nwm.nwm_id(ind);
        k = k+1;
    else
        omit_ind = [omit_ind,i];
    end
end
data(omit_ind,:) = [];
data.nwm_id = nwm_id';

maxes = max(table2array(data(:,19:32)),[],2);
means = nanmean(table2array(data(:,19:32)),2);

subplot(2,1,1)
hist(maxes(means>-500))
xlabel('Maximum nRMSE')
ylabel('Count')

subplot(2,1,2)
hist(means(means>-500))
xlabel('Mean nRMSE')
ylabel('Count')