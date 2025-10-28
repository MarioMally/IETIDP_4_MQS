classdef FullSystemDD < handle
    
    % Class to assemble and solve the DD system of the form:
    % --                -- --    --   --           --
    % |Kpp Kpr  0   Bpp' | |  ap  |   | fp - Kpe*ge |
    % |Krp Krr Brr'  0   | |  ar  | = | fr - Kre*ge |
    % | 0  Brr  0    0   | | lamr |   | gr - Bre*ge |
    % |Bpp  0   0    0   | | lamp |   | gp - Bpe*ge |
    % --                -- --    --   --           --
    % with ae being eliminated/prescribed DOFs ae = ge.
    % The dual-primal approach (based on nullspace method) is
    % taken for solving the system, i.e., we solve
    % --                               -- --    --   --                                  --
    % |Cnull'*Kpp*Cnull Cnull'*Kpr  0   | |  p   |   | Cnull'*(fp - Kpe*ae - Kpp*Corth*q) |
    % |       Krp*Cnull     Krr    Brr' | |  ar  | = |      fr - Kre*ae - Krp*Corth*q     |
    % |        0            Brr     0   | | lamr |   |            gr - Bre*ae             |
    % --                               -- --    --   --                                  --
    % with ap = Cnull*p + Corth*q and (Bpp*Corth)*q = gp - Bpe*ae. Note
    % that ap and ar are split into contributions from ndom subdomains
    % which implies that all matrices and some parts of the RHS can be
    % split into corresponding blocks. These blocks are saved in the
    % property "subdomains" which is of the type Composite<SubSystemDD>.

    properties
        subdomains; % Of type SubSystemDD

        couplRHS; % entails gr and gp
    end

    properties(Access=private)
        nDom; % Number of subsystems
        nMult; % Number of multipliers
        nDOF; % Number of total DOFs
        nPriDOF; % Total number of primal DOFs
        nOrthPriDOF; % Size of q
        nNullPriDOF; % Size of p

        sol; % As cell array of local components
        nullPriSol; % Is p in math notation
        orthPriSol; % Is q in math notation

        multRem; % Remaining multipliers (indices)
        multPri; % Primal multipliers (indices)

        lamRem; % Remaining multipliers
        lamPri; % Primal multipliers

        F; % Matrix representing coarse Problem
        dF; % Factorization of F
        e;
        d;
    end

    methods

        function obj = FullSystemDD(dofSubsets)
            obj.nDom = numel(dofSubsets);
            
            obj.subdomains = cell(obj.nDom,1);
            obj.nDOF = 0;
            obj.nPriDOF = 0;
            for i=1:obj.nDom
                obj.subdomains{i} = SubSystemDD(dofSubsets{i});
                obj.nDOF = obj.nDOF + obj.subdomains{i}.getNumberOfDOFs();
                obj.nPriDOF = obj.nPriDOF + numel(obj.subdomains{i}.getPri());
            end
            
            obj.nMult = obj.subdomains{1}.getNumberOfMults();
            obj.multRem = obj.subdomains{1}.getMultRem();
            obj.multPri = obj.subdomains{1}.getMultPri();

            obj.nOrthPriDOF = numel(obj.multPri);
            obj.nNullPriDOF = obj.nPriDOF - obj.nOrthPriDOF;
        end

        %% Getter functions for outside access of private values
        function res = getSol(obj)
            res = obj.sol;
        end
        function res = getNullPriSol(obj)
            res = obj.nullPriSol;
        end
        function res = getOrthPriSol(obj)
            res = obj.orthPriSol;
        end
        function res = getLamRem(obj)
            res = obj.lamRem;
        end
        function res = getLamPri(obj)
            res = obj.lamPri;
        end
        function res = getNumberOfSubdomains(obj)
            res = obj.nDom;
        end
        function res = getNumberOfDOFs(obj)
            res = obj.nDOF;
        end
        function res = getNumberOfMults(obj)
            res = obj.nMult;
        end
        function res = getSizeOfOrthPriSol(obj)
            res = obj.nOrthPriDOF;
        end
        function res = getSizeOfNullPriSol(obj)
            res = obj.nNullPriDOF;
        end
        function res = getCoarseF(obj)
            res = obj.F;
        end

        %% Setup of subdomains and helper functions
        function setSystemMatrices(obj,KCell)
            for i=1:obj.nDom
                obj.subdomains{i}.setSystemMatrix(KCell{i});
                obj.subdomains{i}.factorizeSystemMatrix();
                obj.subdomains{i}.factorizeForDirichletPrec();
            end           
        end
        function modifySystemMatrices(obj,idxCell,fModCell,modDOFsCell)
            % TBA
            obj.reset();
        end

        function setRHSs(obj,fCell)
            for i=1:obj.nDom
                obj.subdomains{i}.setRHS(fCell{i});
            end
        end
        function modifyRHSs(obj,idxCell,fModCell,modDOFsCell)
            % TBA
            obj.reset();
        end

        function setCouplingMatrices(obj,BCell)
            for i=1:obj.nDom
                obj.subdomains{i}.setCouplingMatrix(BCell{i});
            end
        end
        function modifyCouplingMatrices(obj,idxCell,BModCell,modDOFsCell,modMultsCell)
            % TBA
            obj.reset();
        end

        % function setOrthRepresentations(obj,CorthCell)
        %     for i=1:obj.nDom
        %         obj.subdomains{i}.setOrthRepresentation(CorthCell{i});
        %     end
        % end
        % function modifyOrthRepresentations(obj,idxCell,modCorthCell,modDOFsCell,modPriDOFsCell)
        %     % TBA
        %     obj.reset();
        % end

        function setNullRepresentations(obj,CnullCell)
            for i=1:obj.nDom
                obj.subdomains{i}.setNullRepresentation(CnullCell{i});
            end
        end
        function modifyNullRepresentations(obj,idxCell,modCnullCell,modDOFsCell,modPriDOFsCell)
            % TBA
            obj.reset();
        end

        function setCouplingRHS(obj,g)
            obj.couplRHS = g;
        end
        function modifyCouplingRHS(obj,gMod,modMults)
            % TBA
            obj.reset();
        end

        function setEliminatedRHSs(obj,gEliCell)
            for i=1:obj.nDom
                obj.subdomains{i}.setEliRHS(gEliCell{i});
            end
        end
        function modifyEliminatedRHSs(obj,idxCell,gEliModCell,modDOFsCell)
            % TBA
            obj.reset();
        end

        %% Utilities for the coarse problem
        function multiplyCoarseF(obj)
            if ~isempty(obj.F)
                error('Matrix F is already assembled. Please call obj.reset() to be able to reassamble!');
            end

            obj.F = sparse(obj.nNullPriDOF,obj.nNullPriDOF);
            Id = speye(obj.nNullPriDOF,obj.nNullPriDOF);
            for i=1:obj.nDom
                obj.F = obj.F - obj.subdomains{i}.multiplyF(Id);
            end
        end
        function assembleCoarseF(obj)
            if ~isempty(obj.F)
                error('Matrix F is already assembled. Please call obj.reset() to be able to reassamble!');
            end

            obj.F = sparse(obj.nNullPriDOF,obj.nNullPriDOF);
            Id = speye(obj.nNullPriDOF,obj.nNullPriDOF);
            for i=1:obj.nDom
                obj.F = obj.F - obj.subdomains{i}.multiplyF(Id);
            end
        end
        function factorizeCoarseF(obj)
            if isempty(obj.F)
                error('It is necessary to first assemble F before computing its factorization!');
            end
            obj.dF = decomposition(obj.F,'chol','upper');
        end

        %% Solver utilities
        function computeOrthPriSol(obj)
            % TODO: Think about replacing with iterative method later
            Borth = sparse(obj.nOrthPriDOF,obj.nOrthPriDOF);
            Id = speye(obj.nOrthPriDOF,obj.nOrthPriDOF);
            rhs = obj.couplRHS(obj.multPri);
            for i=1:obj.nDom
                Borth = Borth + obj.subdomains{i}.multiplyOrthPriCoupling(Id);
                rhs = rhs + obj.subdomains{i}.getModificationCouplRHS(obj.multPri);
            end
            obj.orthPriSol = Borth\rhs;
        end

        function computeE(obj)
            if isempty(obj.orthPriSol)
                error('It is necessary to compute orthPriSol before computing e!');
            end

            grModCell = cellfun(@(sub) sub.getModificationCouplRHS(obj.multRem), obj.subdomains, 'UniformOutput', false);
            gr = obj.couplRHS(obj.multRem) + sum([grModCell{:}],2);

            eCell = cellfun(@(sub) sub.computeE(obj.orthPriSol,gr), obj.subdomains, 'UniformOutput', false);
            obj.e = sum([eCell{:}],2);
        end

        function computeD(obj)
            if isempty(obj.orthPriSol)
                error('It is necessary to compute orthPriSol before computing d!');
            end
            dCell = cellfun(@(sub) sub.computeD(obj.orthPriSol), obj.subdomains, 'UniformOutput', false);
            obj.d = sum([dCell{:}],2);
        end

        function res = multiplyInvF(obj,x)
            res = obj.dF\x;
            % [res,~,~,iter] = pcg_w_eigest(-obj.F,x,1e-12,1e4);
        end
        
        function [eigEst,iter] = computeLamRem(obj,tol,maxIter)
            if isempty(obj.e)
                error('It is necessary to compute e before recovering the rem. multipliers!');
            elseif isempty(obj.d)
                error('It is necessary to compute d before recovering the rem. multipliers!');
            elseif isempty(obj.dF)
                error('It is necessary to factorize F first!');
            end
            
            v = obj.multiplyInvF(obj.d);
            GTvCell = cellfun(@(sub) sub.multiplyG(v), obj.subdomains, 'UniformOutput', false);
            rhs = obj.e + sum([GTvCell{:}],2);
            [obj.lamRem, ~, ~, iter, ~, eigEst] = pcg_w_eigest(@obj.evaluateInterfaceProblem,rhs,tol,maxIter,@obj.multiplyDirPrec);
        end

        function res = evaluateInterfaceProblem(obj,x)
            HxCell = cellfun(@(sub) sub.multiplyH(x), obj.subdomains, 'UniformOutput', false);
            res = sum([HxCell{:}],2);
            
            if obj.nOrthPriDOF>0
                GTxCell = cellfun(@(sub) sub.multiplyGT(x), obj.subdomains, 'UniformOutput', false);
                GTx = sum([GTxCell{:}],2);
                v = obj.multiplyInvF(GTx);
                GvCell = cellfun(@(sub) sub.multiplyG(v), obj.subdomains, 'UniformOutput', false);
                
                res = res + sum([GvCell{:}],2);
            end
        end
        
        function recoverNullPriSol(obj)
            if isempty(obj.lamRem)
                error('It is necessary to compute the rem. multipliers before recovering nullPriSol!');
            end
            GTxCell = cellfun(@(sub) sub.multiplyGT(obj.lamRem), obj.subdomains, 'UniformOutput', false);
            obj.nullPriSol = -obj.multiplyInvF(obj.d - sum([GTxCell{:}],2));
        end

        function recoverFullSol(obj)
            if isempty(obj.nullPriSol)
                error('It is necessary to compute nullPriSol before recovering the full solution!');
            end
            obj.sol = cell(obj.nDom,1);
            for i=1:obj.nDom
                obj.subdomains{i}.recoverSolution(obj.lamRem,obj.nullPriSol,obj.orthPriSol);
                obj.sol{i} = obj.subdomains{i}.getSolution();
            end
        end

        function recoverLamPri(obj)
            if isempty(obj.sol)
                error('It is necessary to compute the full sol. before recovering lamPri!');
            end
            error('For now not implemented!');
        end

        %% Preconditioners
        function precomputeScaleMat(obj,boolVal)
            scaleMat = sparse(numel(obj.multRem),numel(obj.multRem));
            for i=1:obj.nDom
                scaleMat = scaleMat - obj.subdomains{i}.localLumpedScaling();
            end
            for i=1:obj.nDom
                obj.subdomains{i}.precomputeScaledLocalCoupling(scaleMat,boolVal);
            end
        end

        function res = multiplyDirPrec(obj,x)
            PxCell = cellfun(@(sub) sub.multiplyDirPrec(x), obj.subdomains, 'UniformOutput', false);
            res = sum([PxCell{:}],2);
        end
        
        %% Other utilities
        function reset(obj)
            obj.F = [];
            obj.dF = [];
            obj.e = [];
            obj.d = [];

            obj.sol = {};
            obj.nullPriSol = [];
            obj.orthPriSol = [];
            obj.lamRem = [];
            obj.lamPri = [];
        end
    end
end