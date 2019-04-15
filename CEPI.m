function  Result = CEPI(Date,Vars,varargin)
%-------------------------------------------------
% 该函数支持常见极端降水指标的计算
% 输入变量：
%         Date为时间矩阵（日尺度）
%         Vars为日降水变量矩阵（每一列代表一个站点或者格点的日降水序列）
% 输入参数：
%          'Percentage'用于确定极端降水量的百分位阈值，可供选项包括：
%                95，默认选项，即日降水序列的第95百分位值
%                a,日降水序列的第a百分位值，通常∈[90,99.5]
%          'Vrange'用于确定极端降水量的固定阈值，可供选项包括：
%                [0,Inf],默认选项
%                [a,b],确定a≤日降水量＜b的基本单位（1日）为极端降水事件，注意a不应＜0且不应＞b
%          'Wet'用于确定湿润天的最小降水量，可供选项包括： 
%                1,默认选项,即日降水量为1mm
%                a,确定a≤日降水量的基本单位（1日）为湿润天
%          'ContiDays'用于确定统计极端降水的基本单位，可供选项包括：
%                1，默认选项,即1日降水量为基本单位
%                a，确定连续a日的降水量为基本单位
%           'Scale'用于确定统计极端降水的尺度，可供选项包括：
%               'd2y',默认选项,即年尺度
%               'd2m'，月尺度
%           'Mrange'用于确定统计统计极端降水的特定顺序月份
%               []，默认选项,无特定顺序月份
%               [a,b,c...]，确定统计a.b.c等顺序月份的极端降水，参数值∈[1,12]
%                           注意参数的顺序特性，如[12,1,2]将选取前一年的12月与该年
%                           的1、2月份而非该年的1、2、12月份，后者为应写作[1,2,12]
% 输出结果:
%          Result.PRCPTOT为湿润天的[年/月]总降水量（wPtot）和降水日数（wDtot）以及
%                        非湿润天的[年/月]总降水量（dPtot）和降水日数（dDtot）
%          Result.CWD为[年/月]尺度上的最大连续湿润天数（MaxW）和连续湿润天数（ConWD）
%          Result.CDD为[年/月]尺度上的最大连续湿非润天数（MaxD）和连续非湿润天数（ConDD）
%          Result.MaxPofConD为连续n天降水的[年/月]最大值（maxV）和对应日期（maxD）
%          Result.PRCPLEV为固定阈值内的日降水的[年/月]总降水量（Pxx）和总降水日数（Dxx）
%          Result.PRCPEXT超过百分位阈值的日降水的[年/月]总降水量（ExP）和总降水日数（ExD）
%
% 例如：Result = CEPI(Date,Vars,'Vrange',[20,Inf],'Percentage',90,'Wet',1,'ContiDays',5,'Scale','d2y','Mrange',[12,1,2])

%% argcheck
narginchk(0,16);
parNames = {'Mrange','Scale','Wet','Vrange','ContiDays','Percentage'};
defaults = {[],'d2y',1,[0,Inf],1,95};
[vMrange,vScale,vWet,vVrange,vContiDays,vPercentage] = internal.stats.parseArgs(parNames, defaults, varargin{:});
ScaleNames = {'d2m','d2y'};  vScale = internal.stats.getParamVal(vScale,ScaleNames,'''Scale''');
if ~isempty(vMrange)&&any(vMrange>12|vMrange<1),   error(' Mrange no greater than 12 or less than 1');end
if vWet<0 ,error('The value of Wet must be positive '); end
if ~internal.stats.isScalarInt(vPercentage,0,100),  error('Percentage no greater than 100 or less than 0');end
if any(vVrange<0)||vVrange(1)>=vVrange(2)||length(vVrange)~=2,    error('The length of Vrange must be 2 and the first value should be less than the second value and both of them should be larger than 0');end
if ~internal.stats.isScalarInt(vContiDays,1,28),   error('The value of ContiDays must be integer from 1 to 28');end
%% Input select function
oVars = Vars;
if ~isempty(vMrange),[Date,Vars] = Input_select(Date,Vars,vMrange);end
%% total precipitation and number of wet days and dry days 
pDate = Date; dVars = bsxfun(@ge,Vars,vWet);  pVars = Vars.*dVars;
[time,wPtot] = ScaleT(pDate,pVars,vScale,'sum');
[~,wDtot] = ScaleT(pDate,dVars,vScale,'sum');
Result.PRCPTOT.Date = time;    Result.PRCPTOT.wPtot = wPtot;
Result.PRCPTOT.wDtot = wDtot;   
%--------------------------
dVars2 = bsxfun(@lt,Vars,vWet); pVars2 = Vars.*dVars2;
[~,dPtot] = ScaleT(pDate,pVars2,vScale,'sum');
[~,dDtot] = ScaleT(pDate,dVars2,vScale,'sum');
Result.PRCPTOT.dPtot = dPtot;   Result.PRCPTOT.dDtot = dDtot;
%--------------------------
%% maximum precipitation in consecutive appointed days 
pDate = Date; pVars=Vars;
[time,maxV,maxD] = MaxPofConD(pDate,pVars,vContiDays,vScale);
Result.MaxPofConD.Date = time;  Result.MaxPofConD.maxV = maxV;
Result.MaxPofConD.maxD = maxD;    
%% total precipitation and number of extreme precipitation days determined by a fixed threshold
pDate = Date; dVars = zeros(size(Vars)); dVars(vVrange(1)<=Vars&Vars<vVrange(2)) = 1;
[time,Pxx] = ScaleT(pDate,Vars.*dVars,vScale,'sum');
[~,Dxx] = ScaleT(pDate,dVars,vScale,'sum');
Result.PRCPLEV.Date = time;   Result.PRCPLEV.Pxx = Pxx;  Result.PRCPLEV.Dxx = Dxx;
%% total precipitation and number of extreme precipitation days determined by a percentile threshold
pDate = Date;   thr = ThreProcess(oVars,vWet,vPercentage);
dVars = bsxfun(@gt,Vars,reshape(thr,1,length(thr)));
[time,ExP] = ScaleT(pDate,Vars.*dVars ,vScale,'sum'); 
[~,ExD] = ScaleT(pDate,dVars,vScale,'sum');
Result.PRCPEXT.Date = time; Result.PRCPEXT.Thresgold = thr;
Result.PRCPEXT.Pxxp = ExP;  Result.PRCPEXT.Dxxp = ExD;
%% the max number  of consecutive wet and dry days  
pDate = Date; dVars = bsxfun(@ge,Vars,vWet);  
[time,ConMaxWet,ConWet] = ContinuityDays(pDate,dVars,vScale);
Result.CWD.Date = time;   Result.CWD.MaxW = ConMaxWet;
Result.CWD.ConWD = ConWet;
%------------------------------
dVars2 = bsxfun(@lt,Vars,vWet); 
[~,ConMaxDry,ConDdry] = ContinuityDays(pDate,dVars2,vScale);
Result.CDD.Date = time;   Result.CWD.MaxDry = ConMaxDry;
Result.CDD.ConDD = ConDdry;
end

function [outDate,outVars] = Input_select(Date,Vars,vMrange)
if isempty(vMrange)
    outDate = Date;
    outVars = Vars;
else
    ID=[];
    for i=length(vMrange):-1:1
        IX = find(Date(:,2)==vMrange(i));
        if i<length(vMrange)
            IX = IX(IX<max(ID)); ID = ID(ID>min(IX));
            ID = [ID;IX(:)];
        else
            ID =[ID;IX(:)];
        end
    end
    sID = sort(ID,'ascend');
    if isequal(vMrange,sort(vMrange))
        outDate = Date(sID,:);
    else 
        outDate = Date(sID,:);
        ID = find(outDate(:,2)==vMrange(1));
        id1 = [ID(diff([-1;ID])~=1);length(outDate)+1];
        for i=1:length(id1)-1
            outDate(id1(i):id1(i+1)-1,1)=outDate(id1(i),1);
        end
    end
    outVars = Vars(sID,:);
end
end
function [OutDate,OutVars] = ScaleT(Date,Vars,vScale,vMethod)
if strcmpi(vScale,'d2m'),  nScale=2;   
elseif strcmpi(vScale,'d2y'),  nScale=1; 
end
[time,~,class]  = unique(Date(:,1:nScale),'rows');
for i=1:length(time)
    if strcmp(vMethod,'max')||strcmp(vMethod,'min')
    eval(['data(',int2str(i),',:)=',vMethod,'(Vars((class==',int2str(i),'),:),[],1);']);
    else
    eval(['data(',int2str(i),',:)=',vMethod,'(Vars((class==',int2str(i),'),:),1);']);
    end
end
    OutDate = time; OutVars = data; 
end
function [time,maxV,maxD] = MaxPofConD(Date,Vars,vContiDays,vScale)
if strcmpi(vScale,'d2m'),        nScale=2;
elseif strcmpi(vScale,'d2y'),    nScale=1;
end   
[time,~,class] = unique(Date(:,1:nScale),'rows');
for i = 1:length(time)       
    ObjectV = Vars((class==i),:);   
    ObjectD = Date((class==i),:);
    clear ConV;
    for j=1:length(ObjectD)-(vContiDays-1)
         ConV(j,:) = sum(ObjectV(j:j+vContiDays-1,:),1);
    end
    maxV(i,:) = max(ConV,[],1);
    
    for k=1:size(Vars,2)
        maxID = find(ConV(:,k)==maxV(i,k));
        for h=1:length(maxID)
        maxD{i,k}{h} = ObjectD(maxID(h):maxID(h)+vContiDays-1,:);
        end
    end
end
end
function [time,ConMax,ConXD] = ContinuityDays(Date,Vars,vScale)
if strcmpi(vScale,'d2m'),        nScale=2;
elseif strcmpi(vScale,'d2y'),    nScale=1;
end   
[time,~,class]  = unique(Date(:,1:nScale),'rows');
for i = 1:length(time)
    ObjectV = Vars((class==i),:);
    ObjectD = Date((class==i),:);
    for j=1:size(ObjectV,2)
        temp = diff([0;ObjectV(:,j);0]);  ids = find(temp==1);  ide = find(temp==-1); 
        Ds = ide-ids;  conD = unique(Ds);
        for h=1:length(conD)
            ConXD{i,j}{h,1} = conD(h);    ConXD{i,j}{h,2} = sum(Ds==conD(h));
            ConXD{i,j}{h,3} = [ObjectD(ids(Ds==conD(h)),:),ObjectD(ide(Ds==conD(h))-1,:)];
            % special max data
            if h==length(conD)
            ConMax(i,j) = conD(h);
            end
        end
    end
end
end
function threshold = ThreProcess(Vars,vWet,vPercentage)
threshold=zeros(1,size(Vars,2));
 for i=1:size(Vars,2)
     flow = Vars(:,i);
     wetflow = flow((flow>=vWet));
     wetflow = sort(wetflow,'ascend');
     threshold(i) = prctile(wetflow,vPercentage);
    % n = length(wetflow);
     %ID = ceil(n*vPercentage/100); 
     %
     %threshold(i) = wetflow(ID);
end
end

