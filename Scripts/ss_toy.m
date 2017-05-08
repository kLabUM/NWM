clear all
close all
cd('~/Google Drive/Docs Kevin/National Water Model/NWM/Scripts');
load('mar29.mat');
load('dec5.mat');
load('granger.mat');
train = dec5;
test = mar29;
granger_ids = granger;
clearvars mar29 dec5 granger
uniq_id = unique(train.ifis_id);
% load('good.mat');
%%
close all
id_ind = [21,27:30,35,53,55:58,60,68,69,76,78,90,91,93,98,99,100,104,110,113,114,129,133,138,151,162,163,164,183,184,191,192,195,197]';
%%
% i = 22;
% ind = find(mar12.ifis_id==uniq_id(id_ind(i)));
% 
% y = mar12.gage_height(ind);
% y(y<-4) = NaN;
% u = mar12.t0(ind);
% 
% train = iddata(y(1:end-500),u(1:end-500),1);
% train = misdata(train);
% train = detrend(train);
% 
% test = iddata(y(end-499:end),u(end-499:end),1);
% test = misdata(test);
% test = detrend(test);
% 
% disp(uniq_id(id_ind(i)));
% plot(train);
% figure;
% plot(test);
%%
lags_train = zeros(length(uniq_id),1);
covs_train = zeros(length(uniq_id),1);
lags_test = zeros(length(uniq_id),1);
covs_test = zeros(length(uniq_id),1);
granger = zeros(length(uniq_id),1);
bottom_width = zeros(length(uniq_id),1); %collinear with order
elevation = zeros(length(uniq_id),1);
mannings = zeros(length(uniq_id),1);
slope = zeros(length(uniq_id),1);
dynamics = zeros(length(uniq_id),1);
tf21 = zeros(length(uniq_id),1);
tf31 = zeros(length(uniq_id),1);
tf32 = zeros(length(uniq_id),1);
ss1 = zeros(length(uniq_id),1);
ss2 = zeros(length(uniq_id),1);
ss3 = zeros(length(uniq_id),1);
ss4 = zeros(length(uniq_id),1);
oe221 = zeros(length(uniq_id),1);
oe211 = zeros(length(uniq_id),1);
ifis_id = zeros(length(uniq_id),1);
data = table(ifis_id,covs_train,lags_train,covs_test,lags_test,granger,bottom_width,elevation,slope,mannings,dynamics,tf21,tf31,tf32,ss1,ss2,ss3,ss4,oe221,oe211);
granger_ids = categorical(granger_ids.site(granger_ids.value==1));

for i = 1:length(uniq_id)
    ind = find(train.ifis_id==uniq_id(i));
    data.ifis_id(i) = string(uniq_id(i));
    data.bottom_width(i) = train.bottom_width(ind(1));
    data.elevation(i) = train.elevation(ind(1));
    data.mannings(i) = train.mannings(ind(1));
    data.slope(i) = train.slope(ind(1));
    data.dynamics(i) = train.dynamics(ind(1));
%     ind = find(train.ifis_id=='CEDARRV03');


    y = train.gage_height(ind);
    y(y<-2) = NaN;
    y = fillmissing(y,'linear');
    u = train.t0(ind);
    
    ind = find(test.ifis_id==uniq_id(i));
%     ind = find(test.ifis_id=='CEDARRV03');
    y_test = test.gage_height(ind);
    y_test(y_test<-2) = NaN;
    y_test = fillmissing(y_test,'linear');
    u_test = test.t0(ind);

    [c,lg] = xcov(y,u,100,'coeff');
    c = c(101:201); %only consider non-neg lags
    lg = lg(101:201);
    [data.covs_train(i),ind] = max(c); 
    data.lags_train(i) = lg(ind);
    
    [c_test,lg_test] = xcov(y_test,u_test,100,'coeff');
    c_test = c_test(101:201); %only consider non-neg lags
    lg_test = lg_test(101:201);
    [data.covs_test(i),ind] = max(c_test);
    data.lags_test(i) = lg_test(ind);
    
    if ismember(uniq_id(i),granger_ids)
        data.granger(i) = 1;
    end
    
    try
        toy_train = iddata(y,u,1);
        toy_train = misdata(toy_train);
        toy_train = detrend(toy_train);

        toy_test = iddata(y_test,u_test,1);
        toy_test = misdata(toy_test);
        toy_test = detrend(toy_test);
    catch
        data.tf21(i) = NaN;
        data.tf31(i) = NaN;
        data.tf32(i) = NaN;
        data.ss1(i) = NaN;
        data.ss2(i) = NaN;
        data.ss3(i) = NaN;
        data.ss4(i) = NaN;
        data.oe221(i) = NaN;
        data.oe211(i) = NaN;
        continue
    end
    
    % TF
    Options = tfestOptions;                                    
    Options.Display = 'off';                                    
    Options.WeightingFilter = [];                              
    Options.InitialCondition = 'backcast';                     
    try
        np = 2;                                                    
        nz = 1;                                                    
        num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
        den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);                               
        iodValue = 0;                                              
        iodFree = true;                                            
        iodMin = 0;                                                
        iodMax = 30;                                               
        sysinit = idtf(num, den, 0);                               
        iod = sysinit.Structure.ioDelay;                           
        iod.Value = iodValue;                                      
        iod.Free = iodFree;                                        
        iod.Maximum = iodMax;                                      
        iod.Minimum = iodMin;                                      
        sysinit.Structure.ioDelay = iod;                           

        tf21 = tfest(toy_train, sysinit, Options);
        [~,fit,~] = compare(toy_test,tf21);
        data.tf21(i) = fit;
    catch
        data.tf21(i) = NaN;
    end
    
    try
        np = 3;                                                    
        nz = 1;                                                    
        num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
        den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
        sysinit = idtf(num, den, 0);                               
        iod = sysinit.Structure.ioDelay;                           
        iod.Value = iodValue;                                      
        iod.Free = iodFree;                                        
        iod.Maximum = iodMax;                                      
        iod.Minimum = iodMin;                                      
        sysinit.Structure.ioDelay = iod;           
        tf31 = tfest(toy_train, sysinit, Options);
        [~,fit,~] = compare(toy_test,tf31);
        data.tf31(i) = fit;
    catch
        data.tf31(i) = NaN;
    end
    
    try
        np = 3;                                                    
        nz = 2;                                                    
        num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
        den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
        sysinit = idtf(num, den, 0);                               
        iod = sysinit.Structure.ioDelay;                           
        iod.Value = iodValue;                                      
        iod.Free = iodFree;                                        
        iod.Maximum = iodMax;                                      
        iod.Minimum = iodMin;                                      
        sysinit.Structure.ioDelay = iod;           
        tf32 = tfest(toy_train, sysinit, Options);
        [~,fit,~] = compare(toy_test,tf32);
        data.tf32(i) = fit;
    catch
        data.tf32(i) = NaN;
    end
    
    %SS
    try
        ss1 = n4sid(toy_train, 1, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [~,fit,~] = compare(toy_test,ss1);
        data.ss1(i) = fit;
    catch
        data.ss1(i) = NaN;
    end
    
    try
        ss2 = n4sid(toy_train, 2, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [~,fit,~] = compare(toy_test,ss2);
        data.ss2(i) = fit;
    catch
        data.ss2(i) = NaN;
    end
    
    try
        ss3 = n4sid(toy_train, 3, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [~,fit,~] = compare(toy_test,ss3);
        data.ss3(i) = fit;
    catch
        data.ss3(i) = NaN;
    end
    
    try
        ss4 = n4sid(toy_train, 4, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [~,fit,~] = compare(toy_test,ss4);
        data.ss4(i) = fit;
    catch
        data.ss4(i) = NaN;
    end
    
    %OE
    Opt = oeOptions;                     
    Opt.InitialCondition = 'backcast';   
    Opt.Focus = 'simulation';            
    try
        nb = 2;                              
        nf = 2;                              
        nk = 1;                              
        oe221 = oe(toy_train,[nb nf nk], Opt);
        [~,fit,~] = compare(toy_test,oe221);
        data.oe221(i) = fit;
    catch
        data.oe221(i) = NaN;
    end
    
    try
        nb = 2;                              
        nf = 1;                              
        nk = 1;                              
        oe211 = oe(toy_train,[nb nf nk], Opt);
        [~,fit,~] = compare(toy_test,oe211);
        data.oe211(i) = fit;
    catch
        data.oe211(i) = NaN;
    end
    
    fprintf('%i/%i\r\n',i,length(uniq_id))
    
    
    
    subplot(2,2,1)
    plot(lg,c)
    
    subplot(2,2,3)
    yyaxis left
    plot(u)
    yyaxis right
    plot(y)
    legend('Flow','Height')
    
    subplot(2,2,2)
    plot(lg_test,c_test)
    
    subplot(2,2,4)
    yyaxis left
    plot(u_test)
    yyaxis right
    plot(y_test)
    legend('Flow','Height')
    
    suptitle(sprintf('%s, Granger=%i',uniq_id(i),data.granger(i)))
    
    pause
end

%%

for i=1:length(good)
% i=81;
    ind = find(train.ifis_id==good(i));
    y_train = train.gage_height(ind);
    y_train(y_train<-2) = NaN;
    %y_train = fillmissing(y_train,'linear');
    u_train = train.t0(ind);
    
    toy_train = iddata(y_train,u_train,1);
    toy_train = misdata(toy_train);
    toy_train = detrend(toy_train);
    
    ind = find(test.ifis_id==good(i));
    y_test = test.gage_height(ind);
    y_test(y_test<-2) = NaN;
    %y_test = fillmissing(y_test,'linear');
    u_test = test.t0(ind);
    
    toy_test = iddata(y_test,u_test,1);
    toy_test = misdata(toy_test);
    toy_test = detrend(toy_test);
    
    % TF
    Options = tfestOptions;                                    
    Options.Display = 'off';                                    
    Options.WeightingFilter = [];                              
    Options.InitialCondition = 'backcast';                     

    np = 2;                                                    
    nz = 1;                                                    
    num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
    den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);                               
    iodValue = 0;                                              
    iodFree = true;                                            
    iodMin = 0;                                                
    iodMax = 30;                                               
    sysinit = idtf(num, den, 0);                               
    iod = sysinit.Structure.ioDelay;                           
    iod.Value = iodValue;                                      
    iod.Free = iodFree;                                        
    iod.Maximum = iodMax;                                      
    iod.Minimum = iodMin;                                      
    sysinit.Structure.ioDelay = iod;                           
           
    tf21 = tfest(toy_train, sysinit, Options);
    
    np = 3;                                                    
    nz = 1;                                                    
    num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
    den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
    sysinit = idtf(num, den, 0);                               
    iod = sysinit.Structure.ioDelay;                           
    iod.Value = iodValue;                                      
    iod.Free = iodFree;                                        
    iod.Maximum = iodMax;                                      
    iod.Minimum = iodMin;                                      
    sysinit.Structure.ioDelay = iod;           
    tf31 = tfest(toy_train, sysinit, Options);
    
    np = 3;                                                    
    nz = 2;                                                    
    num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
    den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
    sysinit = idtf(num, den, 0);                               
    iod = sysinit.Structure.ioDelay;                           
    iod.Value = iodValue;                                      
    iod.Free = iodFree;                                        
    iod.Maximum = iodMax;                                      
    iod.Minimum = iodMin;                                      
    sysinit.Structure.ioDelay = iod;           
    tf32 = tfest(toy_train, sysinit, Options);
    
%     np = 4;                                                    
%     nz = 1;                                                    
%     num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
%     den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
%     sysinit = idtf(num, den, 0);                               
%     iod = sysinit.Structure.ioDelay;                           
%     iod.Value = iodValue;                                      
%     iod.Free = iodFree;                                        
%     iod.Maximum = iodMax;                                      
%     iod.Minimum = iodMin;                                      
%     sysinit.Structure.ioDelay = iod;           
%     tf41 = tfest(toy_train, sysinit, Options);
%     
%     np = 4;                                                    
%     nz = 2;                                                    
%     num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
%     den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
%     sysinit = idtf(num, den, 0);                               
%     iod = sysinit.Structure.ioDelay;                           
%     iod.Value = iodValue;                                      
%     iod.Free = iodFree;                                        
%     iod.Maximum = iodMax;                                      
%     iod.Minimum = iodMin;                                      
%     sysinit.Structure.ioDelay = iod;           
%     tf42 = tfest(toy_train, sysinit, Options);
%     
%     np = 4;                                                    
%     nz = 3;                                                    
%     num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
%     den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
%     sysinit = idtf(num, den, 0);                               
%     iod = sysinit.Structure.ioDelay;                           
%     iod.Value = iodValue;                                      
%     iod.Free = iodFree;                                        
%     iod.Maximum = iodMax;                                      
%     iod.Minimum = iodMin;                                      
%     sysinit.Structure.ioDelay = iod;           
%     tf43 = tfest(toy_train, sysinit, Options);
%     
%     np = 5;                                                    
%     nz = 1;                                                    
%     num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
%     den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
%     sysinit = idtf(num, den, 0);                               
%     iod = sysinit.Structure.ioDelay;                           
%     iod.Value = iodValue;                                      
%     iod.Free = iodFree;                                        
%     iod.Maximum = iodMax;                                      
%     iod.Minimum = iodMin;                                      
%     sysinit.Structure.ioDelay = iod;           
%     tf51 = tfest(toy_train, sysinit, Options);
%     
%     np = 5;                                                    
%     nz = 2;                                                    
%     num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
%     den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
%     sysinit = idtf(num, den, 0);                               
%     iod = sysinit.Structure.ioDelay;                           
%     iod.Value = iodValue;                                      
%     iod.Free = iodFree;                                        
%     iod.Maximum = iodMax;                                      
%     iod.Minimum = iodMin;                                      
%     sysinit.Structure.ioDelay = iod;           
%     tf52 = tfest(toy_train, sysinit, Options);
%     
%     np = 5;                                                    
%     nz = 3;                                                    
%     num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
%     den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
%     sysinit = idtf(num, den, 0);                               
%     iod = sysinit.Structure.ioDelay;                           
%     iod.Value = iodValue;                                      
%     iod.Free = iodFree;                                        
%     iod.Maximum = iodMax;                                      
%     iod.Minimum = iodMin;                                      
%     sysinit.Structure.ioDelay = iod;           
%     tf53 = tfest(toy_train, sysinit, Options);
%     
%     np = 5;                                                    
%     nz = 4;                                                    
%     num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
%     den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
%     sysinit = idtf(num, den, 0);                               
%     iod = sysinit.Structure.ioDelay;                           
%     iod.Value = iodValue;                                      
%     iod.Free = iodFree;                                        
%     iod.Maximum = iodMax;                                      
%     iod.Minimum = iodMin;                                      
%     sysinit.Structure.ioDelay = iod;           
%     tf54 = tfest(toy_train, sysinit, Options);
    
    %SS                                                
    ss1 = n4sid(toy_train, 1, 'DisturbanceModel', 'none', 'Ts', 0, Options);
    ss2 = n4sid(toy_train, 2, 'DisturbanceModel', 'none', 'Ts', 0, Options);
    ss3 = n4sid(toy_train, 3, 'DisturbanceModel', 'none', 'Ts', 0, Options);
    ss4 = n4sid(toy_train, 4, 'DisturbanceModel', 'none', 'Ts', 0, Options);
    
    %OE
    Opt = oeOptions;                     
    Opt.InitialCondition = 'backcast';   
    Opt.Focus = 'simulation';            
    nb = 2;                              
    nf = 2;                              
    nk = 1;                              
    oe221 = oe(toy_train,[nb nf nk], Opt);
    
    nb = 2;                              
    nf = 1;                              
    nk = 1;                              
    oe211 = oe(toy_train,[nb nf nk], Opt);
    
    [~,fit,~] = compare(toy_test,tf21,tf31,tf32,ss1,ss2,ss3,ss4,oe221,oe211);
    fit = fit';
%     compare(toy_test,tf1,tf2,tf3,ss1,ss2,ss3,ss4,oe221,oe211);
    disp(good(i))
    FIT(i,:) = fit;
    fprintf('%i/%i\r\n',i,length(good))
end

%%
for i = 1:10
    ttest(FIT(:,6),FIT(:,i))
end
%%
figure(1)
[~,maxes] = max(FIT');
hist(maxes)
title('Histogram of best model')

figure(2)
[~,mins] = min(FIT');
hist(mins)
title('Histogram of worst model')

figure(3)
boxplot(FIT)
ylim([-50,100])

%%
for i=11:length(good)
% i=81;
    ind = find(train.ifis_id==good(i));
    y_train = train.gage_height(ind);
    y_train(y_train<-2) = NaN;
    %y_train = fillmissing(y_train,'linear');
    u_train = train.t0(ind);
    
    toy_train = iddata(y_train,u_train,1);
    toy_train = misdata(toy_train);
    toy_train = detrend(toy_train);
    
    ind = find(test.ifis_id==good(i));
    y_test = test.gage_height(ind);
    y_test(y_test<-2) = NaN;
    %y_test = fillmissing(y_test,'linear');
    u_test = test.t0(ind);
    
    toy_test = iddata(y_test,u_test,1);
    toy_test = misdata(toy_test);
    toy_test = detrend(toy_test);
    
    % TF
    Options = tfestOptions;                                    
    Options.Display = 'off';                                    
    Options.WeightingFilter = [];                              
    Options.InitialCondition = 'backcast';                               
    iodValue = 0;                                              
    iodFree = true;                                            
    iodMin = 0;                                                
    iodMax = 30;                                               

    np = 4;                                                    
    nz = 3;                                                    
    num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
    den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
    sysinit = idtf(num, den, 0);                               
    iod = sysinit.Structure.ioDelay;                           
    iod.Value = iodValue;                                      
    iod.Free = iodFree;                                        
    iod.Maximum = iodMax;                                      
    iod.Minimum = iodMin;                                      
    sysinit.Structure.ioDelay = iod;           
    tf = tfest(toy_train, sysinit, Options);
    
    
    
    compare(toy_test,tf);
    
    fprintf('%i/%i\r\n',i,length(good))
    pause
end