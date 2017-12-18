clear all
close all
cd('~/Google Drive/Docs Kevin/National Water Model/NWM/Scripts');
load('mar29.mat');
% load('may09.mat');
load('dec5.mat');
load('granger.mat');
train = dec5;
test = mar29;
% test = Feb25toMay09;
granger_ids = granger;
clearvars mar29 dec5 granger Feb25toMay09
uniq_id = unique(train.ifis_id);
data = readtable('IFIS_NWM_data_new.txt');
diff = setdiff(data.ifis_id,uniq_id);
ind = find(data.ifis_id~=diff);
data = data(ind,:);
data.rating = zeros(length(uniq_id),1);

for i = 1:length(data.ifis_id)
    ind = find(train.ifis_id==data.ifis_id(i));
    
    y = train.gage_height(ind);
    u = train.t0(ind);
    
    
    
    ind = find(test.ifis_id==data.ifis_id(i));
    y_test = test.gage_height(ind);
    u_test = test.t0(ind);
    
    y(y<-2) = NaN;
    y_test(y_test<-2) = NaN;
    mini = min([y;y_test]);
    
    y = y - mini+.01;
    y = fillmissing(y,'linear');
    y_test = y_test - mini+.01;
    y_test = fillmissing(y_test,'linear');
    

    try
        model = fit(log(u),log(y),'poly1');
        yhat = exp(feval(model,log(u_test)));
        data.rating(i) = goodnessOfFit(yhat,y_test,'NRMSE');
    catch
        data.rating(i) = NaN;
        continue
    end
end

%%
rating_data = data;

load ss_data_TF_only.mat

data.rating = rating_data.rating;

%%
maxes = nanmax(table2array(data(:,19:32)),[],2);

subplot(1,2,1)
scatter(data.mean,data.rating*100);
axis([-100 100 -100 100])
title('Average nRMSE vs Rating Curve nRMSE')
xlabel('Average nRMSE')
ylabel('Rating Curve nRMSE')
refline(1,0)

subplot(1,2,2)
scatter(maxes,data.rating*100);
axis([-100 100 -100 100])
title('Max nRMSE vs Rating Curve nRMSE')
xlabel('Max nRMSE')
ylabel('Rating Curve nRMSE')
refline(1,0)

%%
for i = 1:9
    subplot(3,3,i)
    temp = table2array(data(:,i+18));
    temp = temp(temp>-100);
    hist(temp)
    title(data.Properties.VariableNames(18+i))
end

%%
bins = [-100:15:100];
[N,X] = hist([maxes,data.mean,data.rating*100],bins);
bar(X,N,1,'grouped')
colormap gray
xlim([-110,100])
legend('Best Transfer Function Model','Average Transfer Function Model','Rating Curve Model')
xlabel('nRMSE')
ylabel('Count')
set(gca,'FontSize',20)
%%
hold on
histogram(maxes,'BinWidth',15)
histogram(data.mean,'BinWidth',15)
histogram(data.rating*100,'BinWidth',15)
legend('Max nRMSE from TF','Max nRMSE from TF','Rating Curve nRMSE','Location','NW')
xlabel('nRMSE')
ylabel('Count')
xlim([-100,100])
hold off