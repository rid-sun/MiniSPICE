#include "minispice.hpp"

bool convNR(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis);


#ifdef MINISPICEKLU
// ref: ../3rd_lib/SuiteSparse-7.7.0/KLU/Tcov/klutest.c
int solve(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int maxIter, int startIdx, int myid) {

    klu_symbolic *&Symbolic = matrix->kluMatrix->KLUmatrixSymbolic;
    klu_numeric *&Numeric = matrix->kluMatrix->KLUmatrixNumeric;
    klu_common &Common = matrix->kluMatrix->KLUmatrixCommon;
    
    int n = matrix->n;
    int *Ap = matrix->kluMatrix->KLUmatrixAp;
    int *Ai = matrix->kluMatrix->KLUmatrixAi;
    double *Ax, *b;

    int ok, save;

    // 由于相似的稀疏性，只是数值不同，所以只做一次分析
    if (startIdx == 0) {
        klu_defaults(&Common);
        Symbolic = klu_analyze(n, Ap, Ai, &Common);
        if (Symbolic == NULL) {
            std::cerr << "Symbolic is null!" << std::endl;
            return MINISPICEERROR;
        }
    }

    int curIter = 0;
    while (++curIter <= maxIter) {

        // laod
        int error = load(matrix, netlist, analysis);
        if (error == MINISPICEERROR) return MINISPICEERROR;

        Ax = matrix->kluMatrix->KLUmatrixAx;
        b = matrix->RHS;

        // first factor
        if (startIdx == 0 && curIter == 1) {
            Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common) ;
        }

        // refactor
        save = Common.scale;
        if (Common.scale == 0) Common.scale = 1;
        klu_refactor(Ap, Ai, Ax, Symbolic, Numeric, &Common);
        Common.scale = save;

        if (Common.status == KLU_SINGULAR) {
            std::cerr << "# singular column : \"" << Common.singular_col << "\"" << std::endl;
        }

        // diagnostics
        klu_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, &Common);
        klu_condest(Ap, Ax, Symbolic, Numeric, &Common);
        klu_rcond(Symbolic, Numeric, &Common);
        klu_flops(Symbolic, Numeric, &Common);

        //
        if (Numeric == NULL || Common.status < KLU_OK) {
            std::cerr << "Common->status is not ok!" <<std::endl;
            return MINISPICEERROR; 
        }

        // solve
        klu_solve(Symbolic, Numeric, n, 1, b, &Common);

        // X_k - X_k_1 = RHS, X_k = oldRHS, X_k_1 = oldRHS - RHS
        // NodeHead &nodeHead = netlist->getNodeHead();
        for (int i = 0; i < matrix->n; i++) {
            matrix->RHS[i] = matrix->oldRHS[i] - matrix->RHS[i];
            // std::cout << "X[" << nodeHead.getID2Name(i + 1) << "]=" << matrix->RHS[i] << " ";
        }
        
        //
        b = matrix->oldRHS;
        matrix->oldRHS = matrix->RHS;
        matrix->RHS = b;

        // 判断收敛性
        if (convNR(matrix, analysis)) break;
    }
    
    return curIter;
}

#endif

#ifdef MINISPICESUPERLU
// 需要显示地初始化mpi环境
// 需要提前进行 superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
// 需要在外层执行
// superlu_gridexit(&grid);
// MPI_Finalize();
int solve(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int maxIter, int startIdx, int myid) {

    superlu_dist_options_t &options = matrix->superluMatrix->options;
    SuperLUStat_t &stat = matrix->superluMatrix->stat;
    SuperMatrix &A = matrix->superluMatrix->A;
    dScalePermstruct_t &ScalePermstruct = matrix->superluMatrix->ScalePermstruct;
    dLUstruct_t &LUstruct = matrix->superluMatrix->LUstruct;
    dSOLVEstruct_t &SOLVEstruct = matrix->superluMatrix->SOLVEstruct;
    gridinfo_t &grid = matrix->superluMatrix->grid;
    int info;
    
    // 

    // 退出标记
    int local_exit_flag = 0;
    int global_exit_flag = 0;

    int curIter = 0;
    while (++curIter <= maxIter) {
        if (grid.iam == 0) {
            // laod
            int error = load(matrix, netlist, analysis);
            if (error == MINISPICEERROR) {
                std::cerr << "Process " << grid.iam << " encountered an error!" << std::endl;
                std::cerr << "load matrix failed!" << std::endl;
                // MPI_Abort(grid->comm, 1);
                local_exit_flag = 1;
            }
            matrix->superluMatrix->global_rhs = matrix->RHS;
        }
        //
        MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, grid.comm);
        if (global_exit_flag) {
            return MINISPICEERROR;
        }

        dcreate_matrix(matrix->superluMatrix, startIdx + curIter);

        if (startIdx + curIter == 1) {
            set_default_options_dist(&options);
            options.IterRefine = NOREFINE;
            options.PrintStat = NO;
            dScalePermstructInit(A.nrow, A.ncol, &ScalePermstruct);
            dLUstructInit(A.ncol, &LUstruct);
            PStatInit(&stat);
        } else {
            options.Fact = SamePattern;
            PStatClear(&stat);
        }

        pdgssvx(&options, &A, &ScalePermstruct, matrix->superluMatrix->local_rhs, matrix->superluMatrix->local_m, 1, &grid,
            &LUstruct, &SOLVEstruct, matrix->superluMatrix->berr, &stat, &info);
        
        if (info) {  /* Something is wrong */
            if (grid.iam == 0) {
                std::cerr << "MINISPICEERROR: INFO = " << info << " returned from pdgssvx()" << std::endl;
                MPI_Abort(grid.comm, 1);
            }
        }
        

        dDestroy_LU(A.ncol, &grid, &LUstruct);

        // 解回收
        MPI_Gatherv(matrix->superluMatrix->local_rhs, matrix->superluMatrix->local_m, MPI_DOUBLE,
                matrix->superluMatrix->global_rhs, matrix->superluMatrix->sendcounts_csr_rhs, matrix->superluMatrix->displs_csr_rhs, MPI_DOUBLE, 0, grid.comm);

        //
        if (grid.iam == 0) {
            // NodeHead &nodeHead = netlist->getNodeHead();
            for (int i = 0; i < matrix->n; i++) {
                matrix->RHS[i] = matrix->oldRHS[i] - matrix->RHS[i];
                // std::cout << "X_" << nodeHead.getID2Name(i + 1) << "=" << matrix->RHS[i] << " ";
            }
            // std::cout << std::endl;
            //
            matrix->superluMatrix->global_rhs = matrix->oldRHS;
            matrix->oldRHS = matrix->RHS;
            matrix->RHS = matrix->superluMatrix->global_rhs;

            // 判断收敛性
            if (convNR(matrix, analysis)) {
                local_exit_flag = 1;
            }
        }

        //
        MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, grid.comm);

        if (global_exit_flag) {
            break;
        }

    }
    
    return curIter;
}

#endif

#ifdef MINISPICEMUMPS

int solve(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int maxIter, int startIdx, int myid) {

    DMUMPS_STRUC_C &id = matrix->mumpsMatrix->id;

    // 退出标记
    int local_exit_flag = 0;
    int global_exit_flag = 0;

    if (startIdx == 0) {
        id.job = JOB_INIT;
        id.par = 1;
        id.sym = 0;
        id.comm_fortran = USE_COMM_WORLD;
        dmumps_c(&id);
        if (id.infog[0] < 0) {
            std::cerr << " (PROC " << myid << " MINISPICEERROR RETURN: \tINFOG(1)= " << id.infog[0] << "\n\t\t\t\t" << "INFOG(2)= " << id.infog[1] << std::endl;
            return MINISPICEERROR;
        }
    }

    int curIter = 0;
    while (++curIter <= maxIter) {

        if (myid == 0) {
            // laod
            int error = load(matrix, netlist, analysis);
            if (error == MINISPICEERROR) {
                std::cerr << "Process " << myid << " encountered an error!" << std::endl;
                std::cerr << "load matrix failed!" << std::endl;
                // MPI_Abort(MPI_COMM_WORLD, 1);
                local_exit_flag = 1;
            }

            id.n = matrix->mumpsMatrix->n;
            id.nnz = matrix->mumpsMatrix->nnz;
            id.irn = matrix->mumpsMatrix->rowidx;
            id.jcn = matrix->mumpsMatrix->colidx;
            id.a = matrix->mumpsMatrix->elementVal;
            id.rhs = matrix->rhs;
        }
        //
        MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (global_exit_flag) {
            return MINISPICEERROR;
        }

        id.ICNTL(1) = -1; 
        id.ICNTL(2) = -1; id.ICNTL(3) = -1; id.ICNTL(4) = 0;
        id.ICNTL(18) = 0;
        id.ICNTL(7) = 5;     //metis
        id.ICNTL(24) = 1;    // 奇异->非奇异

        // 分析
        if (startIdx == 0 && curIter == 1) {
            id.job = 1;
            dmumps_c(&id);
        }

        // LU分解 and 计算
        id.job = 5;
        dmumps_c(&id);

        //
        if (myid == 0) {
            double *tmp = matrix->oldRHS;
            matrix->oldRHS = matrix->RHS;
            matrix->RHS = tmp;

             // 得到新解
            NodeHead &nodeHead = netlist->getNodeHead();  
            for (int i = 0; i < matrix->n; i++) {
                matrix->oldRHS[i] = matrix->RHS[i] - matrix->oldRHS[i];
                std::cout << "X_" << nodeHead.getID2Name(i + 1) << "=" << matrix->oldRHS[i] << " ";
            }
            std::cout << std::endl;

            // 判断收敛性
            if (convNR(matrix, analysis)) local_exit_flag = 1;;
        }

        //
        MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, grid->comm);

        if (global_exit_flag) {
            break;
        }
    }
    
    return curIter;
}

#endif

#ifdef MINSPICEEIGEN

int solve(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int maxIter) {

    
}

#endif

bool convNR(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis) {
    double rtol = analysis->getRtolNR(), atol = analysis->getAtolNR();
    double abs_error = 0.0, rel_error = 0.0;
    for (int i = 0; i < matrix->n; i++) {
        abs_error += (matrix->RHS[i] - matrix->oldRHS[i]) * (matrix->RHS[i] - matrix->oldRHS[i]);
    }
    abs_error = std::sqrt(abs_error);
    // if (abs_error <= atol) return true;
    // else return false;
    for (int i = 0; i < matrix->n; i++) {
        rel_error += matrix->oldRHS[i] * matrix->oldRHS[i];
    }
    rel_error = std::sqrt(rel_error);
    rel_error = abs_error / rel_error;

    if (abs_error <= atol && rel_error <= rtol) return true;
    return false;
}

