function [BTI,Cnull,C_mp2ps] = setupCouplingMatTI(gnum,ndof)
    % mapping from multipatch to patch lvl
    C_mp2ps = patch2multipatch(gnum);

    % Coupling and nullspace for coupling with multiplicity 2
    directCoupled = sum(abs(C_mp2ps),1)==2;
    Cnull1 = C_mp2ps(:,directCoupled);
    [rows1,cols1,data1] = find(Cnull1);
    data1(1:2:end) = -data1(1:2:end);
    B1 = sparse(cols1,rows1,data1,max(cols1),ndof);

    % Coupling and nullspace for coupling with multiplicity >2
    higherCoupling = sum(abs(C_mp2ps),1)>2;
    Cnull2 = C_mp2ps(:,higherCoupling);
    [rows2,cols2,~] = find(Cnull2);

    uniqueCols = unique(cols2);
    
    if ~isempty(uniqueCols)
        rows = [];
        cols = [];
        data = [];
        acc = 0;
        for i=1:numel(uniqueCols)
            col = uniqueCols(i);
            rowsOfFirstEntries = rows2(ismember(cols2,col));
            % Is cicular shift for connecting dofs always correct?
            rowsOfSecondEntries = circshift(rowsOfFirstEntries,1);
            % Remove the last constraint to remove redundancy
            rowsOfFirstEntries = rowsOfFirstEntries(1:end-1);
            rowsOfSecondEntries = rowsOfSecondEntries(1:end-1);
            
            rows = [rows;rowsOfFirstEntries;rowsOfSecondEntries];
            colsOfEntries = acc + (1:numel(rowsOfFirstEntries));
            cols = [cols;colsOfEntries';colsOfEntries'];
            dataOfEntries = ones(numel(rowsOfFirstEntries),1);
            data = [data;-dataOfEntries;dataOfEntries];
    
            acc = acc + numel(rowsOfFirstEntries);
        end
    
        B2 = sparse(cols,rows,data,max(cols),ndof);
    else
        B2 = double.empty(0,size(B1,2));
    end
    
    % Put coupling together
    BTI = [B1;B2];
    Cnull = [Cnull1,Cnull2];

    % Important properties (nullspace representation and full-rank coupling
    assert(nullity(BTI)==0);
    assert(nnz(BTI*Cnull)==0);
end