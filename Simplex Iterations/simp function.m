 function[BFS, A]=simp(A,BV,Cost,Variables)
ZjCj=Cost(BV)*A-Cost;
RUN= true;
while RUN
    ZC=ZjCj(1:end-1);
    if any(ZC<0)
        fprintf('The current BFS is not optimal\n')
        [ent_col,pvt_col]=min(ZC);
        fprintf('Entering Col =%d \n' , pvt_col);
        sol=A(:,end);
        Column=A(:,pvt_col);
        if Column<=0
            error('LPP is unbounded');
        else
            for i=1:size(A,1)
                if Column(i)>0
                    ratio(i)=sol(i)./Column(i);
                else
                    ratio(i)=inf;
                end
            end
            [MinRatio,pvt_row]=min(ratio);
            fprintf('leaving Row=%d \n', pvt_row);
        end
        BV(pvt_row)=pvt_col;
        pvt_key=A(pvt_row,pvt_col);
        A(pvt_row,:)=A(pvt_row,:)./ pvt_key;
        for i=1:size(A,1)
            if i~=pvt_row
                A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
            end
        end
        ZjCj=ZjCj-ZjCj(pvt_col).*A(pvt_row,:)
        ZCj=[ZjCj;A]
        TABLE=array2table(ZCj);
        TABLE.Properties.VariableNames(1:size(ZCj,2))=Variables
        BFS(BV)=A(:,end)
    else
        RUN=false;
        fprintf('Current BFS is Optimal  \n');
        fprintf('Phase End \n')
        BFS=BV;
    end
end
end