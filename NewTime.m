function timeArray = NewTime(T0,T1)
%------------------------------------------------
% 该函数用于生成日期序列，可生成[年/月/日]序列
% 输入变量：
%         T0为日期序列的起始日期
%         T1为日期序列的结束日期
% 输出结果：
%         timeArray为完整的时间序列矩阵
%
% 例如：timeArray = NewTime([1960,2,3],[2018,12,25])
%      生成1960年2月3日至2018年12月25日的完整时间序列矩阵
% v1.0
% copyright@smartcookies

%% argcheck
if ~isequal(size(T0,2),size(T1,2)), error('Please ensure  date ''T0'',''T1'' to same size');end
if ~isequal(fix(T0),T0)||~isequal(fix(T1),T1),error('time nodes must be integer '); 
else
    for i = 1:length(T0)
        switch i
            case 1
                if (T0(i)<0||T1(i)<0),error('Year must larger than 0!');end
            case 2
                if (T0(i)<1||T0(i)>12)||(T1(i)<1||T1(i)>12),error('Month must during 1~12! now !');end
            case 3
                if (T0(i)<1||T0(i)>31)||(T1(i)<1||T1(i)>31),error('Day must during 1~31! now !');end
        end
    end
end
%% main function
Time = [T0;T1];
nDate = size(Time,2);
switch nDate
case 3
IX=0;
while Time(1,1)<=Time(2,1)
    switch Time(1,2)
        case {1 3 5 7 8 10 12}
            d=31;
        case {4 6 9 11}
            d=30;
        case {2}
            if mod(Time(1,1),100)==0&&mod(Time(1,1),400)==0
               d=29;
            elseif mod(Time(1,1),100)~=0&&mod(Time(1,1),4)==0
               d=29;
            else
               d=28;
            end
    end
    timeArray(IX+1:IX+d,3)=1:d;
    timeArray(IX+1:IX+d,2)=Time(1,2);
    timeArray(IX+1:IX+d,1)=Time(1,1);
    IX=IX+d;
    if Time(1,1)==Time(2,1)&&Time(1,2)==Time(2,2)
       break
    end
    Time(1,2)=Time(1,2)+1;
    if Time(1,2)>12
        Time(1,2)=Time(1,2)-12;
        Time(1,1)=Time(1,1)+1;
    end   
end
ID1 = find(timeArray(:,size(Time,2))==Time(1,end),1,'first');
ID2 = find(timeArray(:,size(Time,2))==Time(2,end),1,'last');
timeArray = timeArray(ID1:ID2,:);
case 2
IX=0;
while Time(1,1)<=Time(2,1)
    timeArray(IX+1,2)=Time(1,2);
    timeArray(IX+1,1)=Time(1,1);
    IX=IX+1;
    if Time(1,1)==Time(2,1)&&Time(1,2)==Time(2,2)
       break
    end
    Time(1,2)=Time(1,2)+1;
    if Time(1,2)>12
        Time(1,2)=Time(1,2)-12;
        Time(1,1)=Time(1,1)+1;
    end
 
end          
case 1
timeArray = [Time(1,1):Time(2,1)];
end
end