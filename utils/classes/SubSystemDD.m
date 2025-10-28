classdef SubSystemDD < handle
    properties
        % System of discretized PDE
        K; dKrr; % -> System matrix and its decomposition
        B;
        Cnull; % -> Responsible for coupling
        f; % -> System rhs
        gEli; % -> RHS of elim. constraints (Dirichlet-BC, tree DOFs)
        % Stuff for preconditioning
        dKrVrV;
    end

    % Do that for now because I think the DOF subsets will remain static
    properties(Access=private)

        sol; % -> local solution

        % DOF subsets for TI-DP
        eli;
        rem;
        pri;

        %DOF subsets for Dirichlet precond.
        remVol;
        remInt;

        % Multiplier subsets
        multRem;
        multPri;

        % Some internal variables
        nDOF;
        nMult;

        % Scaling for preconditioning
        scaledB;
    end

    methods
        function obj = SubSystemDD(dofs)
            assert(isempty(intersect(dofs.eli,dofs.pri)),'Eliminated DOFs and primal DOFs are not allowed to intersect!');
            assert(isempty(intersect(dofs.eli,dofs.rem)),'Eliminated DOFs and remaining DOFs are not allowed to intersect!');
            assert(isempty(intersect(dofs.pri,dofs.rem)),'Primal DOFs and remaining DOFs are not allowed to intersect!');

            assert(isempty(intersect(dofs.remVol,dofs.remInt)),'Rem. volume DOFs and rem. interface DOFs are not allowed to intersect!');
            assert(isequal(sort(dofs.rem),union(dofs.remVol,dofs.remInt)),'remVol and remInt are not forming rem together!');

            assert(isempty(intersect(dofs.multPri,dofs.multRem)),'Primal mults and remaining mults are not allowed to intersect!');

            obj.eli = dofs.eli;
            obj.rem = dofs.rem;
            obj.pri = dofs.pri;
            obj.remVol = dofs.remVol;
            obj.remInt = dofs.remInt;

            obj.multRem = dofs.multRem;
            obj.multPri = dofs.multPri;

            obj.nMult = numel(obj.multRem) + numel(obj.multPri);
            obj.nDOF = numel(obj.rem) + numel(obj.pri) + numel(obj.eli);
        end

        %% Getter for DOF and mult subsets
        function eli = getEli(obj)
            eli = obj.eli;
        end
        function rem = getRem(obj)
            rem = obj.rem;
        end
        function pri = getPri(obj)
            pri = obj.pri;
        end
        function remVol = getRemVol(obj)
            remVol = obj.remVol;
        end
        function remInt = getRemInt(obj)
            remInt = obj.remInt;
        end
        function multRem = getMultRem(obj)
            multRem = obj.multRem;
        end
        function multPri = getMultPri(obj)
            multPri = obj.multPri;
        end
        function nDOF = getNumberOfDOFs(obj)
            nDOF = obj.nDOF;
        end
        function nMult = getNumberOfMults(obj)
            nMult = obj.nMult;
        end

        %% Setting of System parts
        function setCouplingMatrix(obj,B)
            if size(B,2)==obj.nDOF && size(B,1)==obj.nMult
                obj.B = B;
            else
                error('Dimension of coupling matrix does not match!');
            end
        end
        function B = getCouplingMatrix(obj)
            B = obj.B;
        end

        function setNullRepresentation(obj,Cnull)
            if size(Cnull,1)==numel(obj.pri)
                obj.Cnull = Cnull;
            else
                error('Size of nullspace representation and number of primal DOFs does not match!');
            end
        end
        function Cnull = getNullRepresentation(obj)
            Cnull = obj.Cnull;
        end

        % function setOrthRepresentation(obj,Corth)
        %     if size(Corth,1)==numel(obj.pri)
        %         obj.Corth = Corth;
        %     else
        %         error('Size of orth. representation and number of primal DOFs does not match!');
        %     end
        % end
        % function Corth = getOrthRepresentation(obj)
        %     Corth = obj.Corth;
        % end

        function setSystemMatrix(obj,K)
            if size(K,1)==size(K,2) && size(K,1)==obj.nDOF
                obj.K = K;
            else
                error('Size of system matrix and number of DOFs does not match!');
            end
        end
        function K = getSystemMatrix(obj)
            K = obj.K;
        end

        function setRHS(obj,f)
            if size(f,1)==obj.nDOF && size(f,2)==1
                obj.f = f;
            else
                error('Size of RHS and number of DOFs does not match!');
            end
        end
        function f = getRHS(obj)
            f = zeros(obj.nDOF,1);
            f(obj.rem) = obj.f(obj.rem) - obj.K(obj.rem,obj.eli)*obj.gEli;
            f(obj.pri) = obj.f(obj.pri) - obj.K(obj.pri,obj.eli)*obj.gEli;
        end

        function setEliRHS(obj,gEli)
            if size(gEli,1)==numel(obj.eli) && size(gEli,2)==1
                obj.gEli = gEli;
            else
                error('Size of RHS and number of DOFs does not match!');
            end
        end
        function gEli = getEliRHS(obj)
            gEli = obj.gEli;
        end

        %% Modify RHS
        function g = getModificationCouplRHS(obj,mults)
            g = - obj.B(mults,obj.eli)*obj.gEli;
        end


        %% Factorizations
        function factorizeSystemMatrix(obj)
            obj.dKrr = decomposition(obj.K(obj.rem,obj.rem),'chol','upper');
        end
        function factorizeForDirichletPrec(obj)
            obj.dKrVrV = decomposition(obj.K(obj.remVol,obj.remVol),'chol','upper');
        end
        
        %% Matrix-Vector products
        % Stiffness Components
        function res = multiplyInvKrr(obj,x)
            res = obj.dKrr\x;
        end
        function res = multiplyKpp(obj,x)
            res = obj.multiplyNullRepT(obj.K(obj.pri,obj.pri)*obj.multiplyNullRep(x));
        end
        function res = multiplyKrp(obj,x)
            res = obj.K(obj.rem,obj.pri)*obj.multiplyNullRep(x);
        end
        function res = multiplyKpr(obj,x)
            res = obj.multiplyNullRepT(obj.K(obj.pri,obj.rem)*x);
        end
        
        % Coupling
        function res = multiplyBr(obj,x)
            res = obj.B(obj.multRem,obj.rem)*x;
        end
        function res = multiplyBrT(obj,x)
            res = obj.B(obj.multRem,obj.rem)'*x;
        end
        function res = multiplyBp(obj,x)
            res = obj.B(obj.multPri,obj.pri)*x;
        end
        function res = multiplyBpT(obj,x)
            res = obj.B(obj.multPri,obj.pri)'*x;
        end

        % Nullspace
        function res = multiplyNullRep(obj,x)
            res = obj.Cnull*x;
        end
        function res = multiplyNullRepT(obj,x)
            res = obj.Cnull'*x;
        end

        % Preconditioning
        function res = multiplyInvKrVrV(obj,x)
            res = obj.dKrVrV\x;
        end
        function res = multiplyKrIrV(obj,x)
            res = obj.K(obj.remInt,obj.remVol)*x;
        end
        function res = multiplyKrVrI(obj,x)
            res = obj.K(obj.remVol,obj.remInt)*x;
        end
        function res = multiplyKrIrI(obj,x)
            res = obj.K(obj.remInt,obj.remInt)*x;
        end


        %% Functionalities related to orth. component (Bp')
        function res = multiplyOrthPriCoupling(obj,x)
            res = obj.multiplyBp(obj.multiplyBpT(x));
        end

        %% Local Schur complement components
        function res = multiplyF(obj,x)
            res = obj.multiplyKpr(obj.multiplyInvKrr(obj.multiplyKrp(x)))...
                  - obj.multiplyKpp(x);
        end

        function res = multiplyG(obj,x)
            res = obj.multiplyBr(obj.multiplyInvKrr(obj.multiplyKrp(x)));
        end

        function res = multiplyGT(obj,x)
            res = obj.multiplyKpr(obj.multiplyInvKrr(obj.multiplyBrT(x)));
        end

        function res = multiplyH(obj,x)
            res = obj.multiplyBr(obj.multiplyInvKrr(obj.multiplyBrT(x)));
        end

        function res = computeD(obj,q)
            frMod = obj.f(obj.rem) - obj.K(obj.rem,obj.eli)*obj.gEli...
                                   - obj.K(obj.rem,obj.pri)*obj.multiplyBpT(q);
            fpMod = obj.f(obj.pri) - obj.K(obj.pri,obj.eli)*obj.gEli...
                                   - obj.K(obj.pri,obj.pri)*obj.multiplyBpT(q);
            res = obj.multiplyKpr(obj.multiplyInvKrr(frMod)) - obj.multiplyNullRepT(fpMod);
        end

        function res = computeE(obj,q,gr)
            frMod = obj.f(obj.rem) - obj.K(obj.rem,obj.eli)*obj.gEli...
                                   - obj.K(obj.rem,obj.pri)*obj.multiplyBpT(q);
            res = obj.multiplyBr(obj.multiplyInvKrr(frMod)) - gr;
        end

        %% Dirichlet preconditioner with super-lumped scaling

        function res = localLumpedScaling(obj)
            A = sparse(1:numel(obj.remInt),1:numel(obj.remInt),1./diag(obj.K(obj.remInt,obj.remInt)));
            res = obj.B(obj.multRem,obj.remInt)*A*obj.B(obj.multRem,obj.remInt)';
        end

        function precomputeScaledLocalCoupling(obj,scaleMat,boolVal)
            if boolVal
                A = sparse(1:numel(obj.remInt),1:numel(obj.remInt),1./diag(obj.K(obj.remInt,obj.remInt)));
                obj.scaledB = scaleMat\(obj.B(obj.multRem,obj.remInt)*A);
            else
                obj.scaledB = obj.B(obj.multRem,obj.remInt);
            end
        end

        function res = multiplyDirPrec(obj, x)
            res = obj.scaledB*(obj.multiplyKrIrI(obj.B(obj.multRem,obj.remInt)'*x)... 
                  - obj.multiplyKrIrV(obj.multiplyInvKrVrV(obj.multiplyKrVrI(obj.scaledB'*x))));
        end

        %% Steps for solution recovery
        function recoverSolution(obj,lam,p,q)
            % This func. assumes that obj.f still contains the modification due
            % to the elimination of DOFs
            obj.sol = zeros(obj.nDOF,1);
            
            obj.sol(obj.eli) = obj.gEli;

            obj.sol(obj.pri) = obj.multiplyNullRep(p) + obj.multiplyBpT(q);
            
            obj.sol(obj.rem) = obj.multiplyInvKrr(...
                obj.f(obj.rem)...
                - obj.K(obj.rem,obj.eli)*obj.gEli...
                - obj.K(obj.rem,obj.pri)*obj.sol(obj.pri)...
                - obj.multiplyBrT(lam)...
            );
        end

        function sol = getSolution(obj)
            sol = obj.sol;
        end

    end
end