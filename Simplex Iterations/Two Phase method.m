format short;
clear all;
clc;
%||| convert to maximization ||
%||| write in standard form |||
%||| introduce dummy variables in equations and subtract them in Z |||

%set this variables after converting equation to standard form
variables = {'x1','x2','x3','s1','s2','A1','A2','sol'}; %variables in Z
ovariables = {'x1','x2','x3','s1','s2','sol'}; %variables without artifical
%set Z values in origC
origC = [-7.5, 3,0,0,0,-1,-1,0];  % -1 are for A1 and A2 0 is for solution
%matrix coefficients in info
info = [3 -1 -1 -1 0 1 0 3;1 -1 1 0 -1 0 1 2]; %matrix for constraints RHS also combined
BV = [6 7] ; %column indices for the identity matrix here A1 and A2

%phase 1 of method
%consider only the cost for A1 and A2 here
cost = [0 0 0 0 0 -1 -1 0];
A = info;
startBV = find(cost<0); %store the artificial variables index as start BV
%compute Zj - Cj
ZjCj = cost(BV)*A - cost;
%put ZjCj row above the A matrix and print as table
initialTable = array2table([ZjCj;A]);
initialTable.Properties.VariableNames(1:size(A,2))=variables

RUN = true;
while RUN
    %find the most negative value excluding the sol column
    ZC = ZjCj(1,1:end-1);
    if any(ZC<0)

        %find entering variable
        [EnterCol pvt_col] = min(ZC);  %most negative value and pivot column index
        fprintf('pivot column is : %d\n',pvt_col);
        %find leaving variable
        %find ratio of last column and pivot column values
        sol = A(:,end); %last column
        column = A(:,pvt_col); %pivot column
        %if all pivot column values are negative then unbounded solution
        if column<0
            fprintf('Unbounded solution\n');
        else
            for i=1:size(A,1) %for every value in pivot column
                if column(i)>0 %if positive then get ratio
                    ratio(i) = sol(i)./column(i);
                else  %else store as inf
                    ratio(i) = inf;
                end
            end
            %find the minimum ratio
            [minRatio,pvt_row] = min(ratio);
            fprintf('pivot row is : %d\n',pvt_row);
        end
        %update basic variables
        BV(pvt_row)  = pvt_col;
        %set pivot element
        pvt_key = A(pvt_row,pvt_col);
        %set pvt_key as 1 and all other entries in column as 0
        A(pvt_row,:) = A(pvt_row,:)./pvt_key; %divide entire row by pvtkey
        for i = 1:size(A,1)   %for all rows other than pivot row
            if i~=pvt_row
                A(i,:) = A(i,:) - A(i,pvt_col).*A(pvt_row,:);  %make element as zero
            end
        end
        %updating the Zj - Cj row
        ZjCj = ZjCj - ZjCj(pvt_col).*A(pvt_row,:);
        %print the table
        table = array2table([ZjCj;A]); 
        table.Properties.VariableNames(1:size(A,2))=variables
        %print the BFS
        BFS(BV) = A(:,end);

    else %if no ZjCj is negative then solution is reached
        RUN = false;
        fprintf('optimal solution is reached for phase one\n');
        BFS = BV;
    end
end

fprintf('PHASE 2 STARTS\n');
%remove artificial variables from the table and cost
A(:,startBV) = [];
origC(:,startBV) = [];
%proceed as normal simplex method
[optBFS, optA] = simp(A,BFS,origC,ovariables);

%printing final solution
final_BFS = zeros(1,size(A,2));
final_BFS(optBFS) = optA(:,end);
final_BFS(end) = sum(final_BFS.*origC);  %computing Z value

finalTable = array2table(final_BFS);
finalTable.Properties.VariableNames(1:size(finalTable,2)) = ovariables




