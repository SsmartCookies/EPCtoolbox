function  Result = CEPI(Date,Vars,varargin)
%-------------------------------------------------
% �ú���֧�ֳ������˽�ˮָ��ļ���
% ���������
%         DateΪʱ������ճ߶ȣ�
%         VarsΪ�ս�ˮ��������ÿһ�д���һ��վ����߸����ս�ˮ���У�
% ���������
%          'Percentage'����ȷ�����˽�ˮ���İٷ�λ��ֵ���ɹ�ѡ�������
%                95��Ĭ��ѡ����ս�ˮ���еĵ�95�ٷ�λֵ
%                a,�ս�ˮ���еĵ�a�ٷ�λֵ��ͨ����[90,99.5]
%          'Vrange'����ȷ�����˽�ˮ���Ĺ̶���ֵ���ɹ�ѡ�������
%                [0,Inf],Ĭ��ѡ��
%                [a,b],ȷ��a���ս�ˮ����b�Ļ�����λ��1�գ�Ϊ���˽�ˮ�¼���ע��a��Ӧ��0�Ҳ�Ӧ��b
%          'Wet'����ȷ��ʪ�������С��ˮ�����ɹ�ѡ������� 
%                1,Ĭ��ѡ��,���ս�ˮ��Ϊ1mm
%                a,ȷ��a���ս�ˮ���Ļ�����λ��1�գ�Ϊʪ����
%          'ContiDays'����ȷ��ͳ�Ƽ��˽�ˮ�Ļ�����λ���ɹ�ѡ�������
%                1��Ĭ��ѡ��,��1�ս�ˮ��Ϊ������λ
%                a��ȷ������a�յĽ�ˮ��Ϊ������λ
%           'Scale'����ȷ��ͳ�Ƽ��˽�ˮ�ĳ߶ȣ��ɹ�ѡ�������
%               'd2y',Ĭ��ѡ��,����߶�
%               'd2m'���³߶�
%           'Mrange'����ȷ��ͳ��ͳ�Ƽ��˽�ˮ���ض�˳���·�
%               []��Ĭ��ѡ��,���ض�˳���·�
%               [a,b,c...]��ȷ��ͳ��a.b.c��˳���·ݵļ��˽�ˮ������ֵ��[1,12]
%                           ע�������˳�����ԣ���[12,1,2]��ѡȡǰһ���12�������
%                           ��1��2�·ݶ��Ǹ����1��2��12�·ݣ�����ΪӦд��[1,2,12]
% ������:
%          Result.PRCPTOTΪʪ�����[��/��]�ܽ�ˮ����wPtot���ͽ�ˮ������wDtot���Լ�
%                        ��ʪ�����[��/��]�ܽ�ˮ����dPtot���ͽ�ˮ������dDtot��
%          Result.CWDΪ[��/��]�߶��ϵ��������ʪ��������MaxW��������ʪ��������ConWD��
%          Result.CDDΪ[��/��]�߶��ϵ��������ʪ����������MaxD����������ʪ��������ConDD��
%          Result.MaxPofConDΪ����n�콵ˮ��[��/��]���ֵ��maxV���Ͷ�Ӧ���ڣ�maxD��
%          Result.PRCPLEVΪ�̶���ֵ�ڵ��ս�ˮ��[��/��]�ܽ�ˮ����Pxx�����ܽ�ˮ������Dxx��
%          Result.PRCPEXT�����ٷ�λ��ֵ���ս�ˮ��[��/��]�ܽ�ˮ����ExP�����ܽ�ˮ������ExD��
%
% ���磺Result = CEPI(Date,Vars,'Vrange',[20,Inf],'Percentage',90,'Wet',1,'ContiDays',5,'Scale','d2y','Mrange',[12,1,2])

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

