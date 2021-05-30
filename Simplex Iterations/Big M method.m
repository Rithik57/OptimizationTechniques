format short;
clear all;
clc;
% ||| convert problem to maximization |||
% ||| write equation in standard form |||
% ||| add artifical variables required to make identity matrix here 2 |||
%%maximize Z = -2x1 - x2 - MA1 - MA2
% 3x1 +x2 + A1 = 3;
% 4x1 + 3x2 -s2 + A2 = 6;
% x1 + 2x2 +s3 = 3;

variables = {'x1','x2','s2','s3','A1','A2','sol'};
M = 1000; %large number
cost = [-2 -1 0  0 -M -M 0]; %last zero for sol
A = [3 1 0 0 1 0 3;  %coefficients matrix
    4 3 -1 0 0 1 6;
    1 2 0 1 0 0 3];
s = eye(size(A,1));
%set initial basic variables as A1 A2 and s3 (according to identity matrix formed)
BV = [5 6 4];


%find Zj - Cj
ZjCj = cost(BV)*A - cost;
Zcj = [ ZjCj ; A ];
simpTable = array2table(Zcj);
simpTable.Properties.VariableNames(1:size( Zcj , 2)) = variables

RUN = true
while RUN

ZC = ZjCj (:,1:end-1); %leave sol column
if any(ZC<0)
    fprintf('current BFS is not optimal\n');
    %find the entering variable and pivot column
    %select min Zj - Cj 
    [entVal , pvt_col] = min(ZC);
    fprintf('pivot column is : %d\n',pvt_col);
    %find pivot row
    %select min ratio
    sol = A(:,end);
    column = A(:,pvt_col);
    if all(column)<=0
        fprintf('solution is unbounded\n')
    else
        for i=1:size(A,1) %for every row
                if column(i)>0   %if pivot column value is pos
                    ratio(i)=sol(i)./column(i);  %find ratio
                else
                    ratio(i)=inf;  %else set to inf
                end
            end
            [minRatio,pvt_row]=min(ratio);  %get the min value
            fprintf('pivot Row is =%d \n', pvt_row); %pvt_row is the index
            %update basic variable and pivot element
            BV(pvt_row) = pvt_col;
            pvt_key = A(pvt_row,pvt_col);
            %update the table for next
            A(pvt_row,:) = A(pvt_row,:)./pvt_key;
            for i=1:size(A,1)
                if i~=pvt_row
                A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
                end
            end
            %calculate Zj - Cj
            ZjCj = cost(BV)*A - cost;
            ZCj = [ZjCj ; A];
            table = array2table(ZCj);
            table.Properties.VariableNames(1:size(ZCj,2)) = variables
    end %if not unbounded
    
else %if optimal
    RUN = false;
    fprintf('BFS is optimal\n');
    BFS = BV;
end
end %while loop ends
%final optimal solution
final_BFS = zeros(1,size(A,2)); %equal to number of columns
final_BFS(BV) = A(:,end);  %update the values of the BV as the ones in the end column
final_BFS(end) = sum(final_BFS.*cost); %find the final value of Z by summing
%and store in the last sol value of BFS
finalTable = array2table(final_BFS);
finalTable.Properties.VariableNames(1:size(A,2)) = variables

