/*
* @Author: SsmartCookies
* @Date:   2019-04-15 17:07:19
* @Last Modified by:   SsmartCookies
* @Last Modified time: 2019-04-15 22:41:58
*/
function [outDate,outVars] = ScaleT(Date,Vars,varargin)
%---------------------------------------------------
% 该函数支持气候变量从高时间分辨率向低时间分辨率转换，最低时间分辨率为日（d）
% 输入变量：
%         Date为时间矩阵（年月日）
%         Vars为气候变量矩阵（支持多维度）
% 输入参数：
%         'Scale'用于确定尺度转换单位，可供选项包括：
%               []，默认选项，不进行尺度转换
%               'd2m'，日尺度 -> 月尺度
%               'd2y'，日尺度 -> 年尺度
%               'm2y'，月尺度 -> 年尺度
%         'Method'用于确定要素转换方法，可供选项包括：
%               'mean'，默认选项，取均值
%               'median'，取中位值
%               'max'，取最大值
%               'min'，取最小值
%               'sum'，取加和值
%         'Mrange'用于确定尺度转换的特定顺序月份，可供选项包括：
%               []，默认选项,无特定顺序月份
%               [a,b,c...]，确定a.b.c等顺序月份的要素进行尺度转换，参数值∈[1,12]
%                           注意参数的顺序特性，如[12,1,2]将选取前一年的12月与该年
%                           的1、2月份而非该年的1、2、12月份，后者为应写作[1,2,12]
%         'Yrange'用于确定尺度转换的特定年份，可供选项包括：
%               []，默认选项,无特定年份
%               [a,b,c...]，确定a.b.c等年份的要素进行尺度转换，超出输入变量（Date）
%                           中的最小/大年份的参数为非法输入参数
%         'Dm'用于确定进行尺度转换的变量所在的维度，可供选项包括：
%               '1'，默认选项,即默认变量值按行排列，变量名按列排列
%               'a',确定第a维为尺度转换的维度，超出输入变量（Vars）的最大纬度的参数
%                   为非法输入参数
% 输出结果：
%         outDate为尺度转换后的时间矩阵
%         outVars为尺度转换后的气候变量矩阵
%
% 例如：[outDate,outVars] = ScaleT(Date,Vars,'Scale','d2y','Method','mean','Mrange',[4:8],'Yrange',[1980:2010],'Dm',1)
%      即：将1980-2010年4月至8月的气候变量数据求取均值作为对应年份的气候特征值
%
% v1.0
% copyright@smartcookies

%% argcheck
narginchk(0,12);
parNames = {'Scale','Method','Dm','Mrange','Yrange'};
default = {[],'mean',1,[],[]};
[vScale,vMethod,vDm,vMrange,vYrange] = internal.stats.parseArgs(parNames,default,varargin{:});
ScaleName = {'d2m','d2y','m2y'};
vScale = internal.stats.getParamVal(vScale,ScaleName,'''Scale''');
MethodNames = {'mean','median','max','min','sum'};
vMethod = internal.stats.getParamVal(vMethod,MethodNames,'''Method''');
if vDm>length(size(Vars))||vDm<0
    error('Dm must greater than 0 and less than %d',length(size(Vars)));end
if ~isempty(vMrange)&&any(vMrange>12|vMrange<1)
   error(' Mrange no greater than 12 or less than 1');end
if ~isempty(vYrange)&&any(vYrange>max(Date(:,1))|vYrange<min(Date(:,1)))
   error(' Yrange must during %d to %d',min(Date(:,1)),max(Date(:,1)));end
if ~isequal(size(Date,1),size(Vars,vDm))
    error('Length of Date and Vars is not equal !');end

%% main function
[Date,Vars] = OrderCheck(Date,Vars);
 if isempty(vMrange)&&isempty(vYrange),   vDate = Date;  vVars = Vars;
 else [vDate,vVars] = Input_select(Date,Vars,vDm,vMrange,vYrange);end
 if isempty(vScale), outDate = vDate; outVars = vVars; return;
 elseif strcmpi(vScale,'d2m')&&isequal(size(vDate,2),3),    nScale = 2;
 elseif strcmpi(vScale,'d2y')&&isequal(size(vDate,2),3),    nScale = 1;
 elseif strcmpi(vScale,'m2y')&&isequal(size(vDate,2),2),    nScale = 1;
 end
[outDate,~,class]  = unique(vDate(:,1:nScale),'rows');
for i=1:length(outDate)
  ID = find(class==i);
  spVars = DimVars(vVars,vDm,ID);
  str = Dmn(Vars,vDm,i);
    if strcmpi(vMethod,'max')||strcmpi(vMethod,'min')
        eval([str,'=',vMethod,'(spVars,[],',int2str(vDm),');']);
    elseif strcmpi(vMethod,'mean')||strcmpi(vMethod,'median')||strcmpi(vMethod,'sum')
        eval([str,'=',vMethod,'(spVars,',int2str(vDm),');']);
    end
end

end
function [Date,Vars] = OrderCheck(Date,Vars)
nDate = size(Date,2);
outDate0 = Date(1,1:nDate);
outDate1 = Date(end,1:nDate);
if nDate==1||nDate==3, Tlong = datenum(outDate1)-datenum(outDate0)+1;
elseif nDate==2,Tlong = sum((datenum(outDate1)-datenum(outDate0)).*[12,1])+1;
end
 ismiss = ~isequal(Tlong,size(Date,1));
if ismiss, warning('There is missing date value in Date !Please check it!');
else
    [~,order,~] = unique(Date,'rows');
    issort = isequal(order,sort(order,'ascend'));
    if issort, return;
    else Date = Date(order,:); Vars = Vars(order,:);
    end
end
end
function [Date,Vars] = Input_select(Date,Vars,vDm,vMrange,vYrange)
if ~isempty(vYrange)
    ID = [];
    for i=1:length(vYrange)
        IX = find(Date(:,1)==vYrange(i));
        ID = [ID;IX(:)];
    end
    sID = sort(ID,'ascend');
    tyDate = Date(sID,:);
    tyVars = DimVars(Vars,vDm,sID);
else
    tyDate = Date;
    tyVars = Vars;
end
if isempty(vMrange)
    Date = tyDate;
    Vars = tyVars;
else
    ID=[];
    for i=length(vMrange):-1:1
        IX = find(tyDate(:,2)==vMrange(i));
        if i<length(vMrange)
            IX = IX(IX<max(ID)); ID = ID(ID>min(IX));
            ID = [ID;IX(:)];
        else
            ID =[ID;IX(:)];
        end
    end
    sID = sort(ID,'ascend');
    if isequal(vMrange,sort(vMrange))
        Date = tyDate(sID,:);
    else
        Date = tyDate(sID,:);
        ID = find(Date(:,2)==vMrange(1));
        id1 = [ID(diff([-1;ID])~=1);length(Date)+1];
        for i=1:length(id1)-1
            Date(id1(i):id1(i+1)-1,1)=Date(id1(i),1);
        end
    end
    Vars = DimVars(tyVars,vDm,sID);
end
end
function spVars = DimVars(Vars,Dm,ID)
ID = reshape(ID,1,length(ID));
nDm = length(size(Vars));
str = 'spVars=Vars(';n = 1;
while n<=nDm
    if n==Dm&&n==nDm
        str = strcat(str,'[',int2str(ID),']',')');
    elseif n==Dm&&n~=nDm
        str = strcat(str,'[',int2str(ID),']',',');
    elseif n==nDm
        str = strcat(str,':)');
    else
        str = strcat(str,':,');
    end
    n = n+1;
end
eval([str,';']);
end
function str = Dmn(Vars,Dm,ni)
ni = reshape(ni,1,length(ni));
nDm = length(size(Vars));
str = 'outVars(';n = 1;
while n<=nDm
    if n==Dm&&n==nDm
        str = strcat(str,int2str(ni),')');
    elseif n==Dm&&n~=nDm
        str = strcat(str,int2str(ni),',');
    elseif n==nDm
        str = strcat(str,':)');
    else
        str = strcat(str,':,');
    end
    n = n+1;
end
end

