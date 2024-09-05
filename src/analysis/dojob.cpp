#include "minispice.hpp"

// 函数前向声明
int homotopy(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int myid);
int newton_raphson(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int maxNRIter, int startIdx, int myid);
int pseudoTran(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int myid);
int tranAn(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int myid);
int ptConv(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis);


int dojob(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, std::string outFileName, int myid) {

    int local_exit_flag = 0;
    int global_exit_flag = 0;

    AnalysisType curmode = analysis->getOpType();
    analysis->setCurMode(curmode);

    // 这里阻塞，共享模式，下面就都是全局的了
    int curmode_value = static_cast<int>(curmode);
   
#ifdef MINISPICESUPERLU
    MPI_Bcast(&curmode_value, 1, MPI_INT, 0, matrix->superluMatrix->grid.comm);
#else
    MPI_Bcast(&curmode_value, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    curmode = static_cast<AnalysisType>(curmode_value);

    analysis->setCurMode(curmode);

    if (curmode == AnalysisType::no_An) return MINISPICEOK;
    /* 当前模式一定显示op */
    if (curmode == AnalysisType::TRAN) return MINISPICEERROR;

MPI_ENV_0_PID_BEGIN(myid)

    switch (analysis->getCurMode()) {
    case AnalysisType::DC_NR:
        std::cout << "---Curmode is DC_NR" << std::endl;
        break;
    case AnalysisType::DC_Homotopy:
        std::cout << "---Curmode is DC_Homotopy" << std::endl;
        break;
    case AnalysisType::DC_PTran:
        std::cout << "---Curmode is DC_Ptran" << std::endl;
        break;
    case AnalysisType::TRAN:
        std::cout << "---Curmode is Tran" << std::endl;
        break;
    case AnalysisType::no_An:
        std::cout << "---Curmode is no_An" << std::endl;
        break;
    default:
        break;
    }

    // 1. setup [创建矩阵、创建映射关系]
    std::cout << "---Setup beginning" << std::endl;
    setup(matrix, netlist, analysis, outFileName);
    std::cout << "---Setup ending" << std::endl;

    switch (analysis->getCurMode()) {
    case AnalysisType::DC_NR:
        std::cout << "---DC_NR analysis beginning" << std::endl;
        break;
    case AnalysisType::DC_Homotopy:
        std::cout << "---DC_Homotopy analysis beginning" << std::endl;
        break;
    case AnalysisType::DC_PTran:
        std::cout << "---DC_Ptran analysis beginning" << std::endl;
        break;
    case AnalysisType::TRAN:
        std::cout << "---DC_Tran analysis beginning" << std::endl;
        break;
    case AnalysisType::no_An:
        std::cout << "---no_An" << std::endl;
        break;
    default:
        break;
    }

MPI_ENV_0_PID_END

    // 2. 执行 DC 分析算法
    int conv = 0;
    switch (curmode) {
    case AnalysisType::DC_Homotopy:
        conv = homotopy(matrix, netlist, analysis, myid);
        break;
    case AnalysisType::DC_NR:
        conv = newton_raphson(matrix, netlist, analysis, analysis->getMaxNRIterNum(), 0, myid);
        break;
    case AnalysisType::DC_PTran:
        conv = pseudoTran(matrix, netlist, analysis, myid);
        break;
    default:
        break;
    }

MPI_ENV_0_PID_BEGIN(myid)

    if (conv == NONCONV) {
        std::cerr << "------DC non-convergence!" << std::endl;
        local_exit_flag = 1;
        // return MINISPICEERROR;
    } else {
        std::cout << "------Total NR_Iters is " << conv << std::endl;
        std::cout << "------DC convergence!" << std::endl;
    }

    switch (analysis->getCurMode()) {
    case AnalysisType::DC_NR:
        std::cout << "---DC_NR analysis ending" << std::endl;
        break;
    case AnalysisType::DC_Homotopy:
        std::cout << "---DC_Homotopy analysis ending" << std::endl;
        break;
    case AnalysisType::DC_PTran:
        std::cout << "---DC_Ptran analysis ending" << std::endl;
        break;
    case AnalysisType::TRAN:
        std::cout << "---DC_Tran analysis ending" << std::endl;
        break;
    case AnalysisType::no_An:
        std::cout << "---no_An" << std::endl;
        break;
    default:
        break;
    }

MPI_ENV_0_PID_END

#ifdef MINISPICESUPERLU
    MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, matrix->superluMatrix->grid.comm);
#else
    MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (global_exit_flag) {
        return MINISPICEERROR;
    }

    unsetup(matrix, netlist, analysis);

    // 3. 执行 Tran 分析算法
    curmode = analysis->getTranType();
    analysis->setCurMode(curmode);

    curmode_value = static_cast<int>(curmode);
#ifdef MINISPICESUPERLU
    MPI_Bcast(&curmode_value, 1, MPI_INT, 0, matrix->superluMatrix->grid.comm);
#else
    MPI_Bcast(&curmode_value, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    curmode = static_cast<AnalysisType>(curmode_value);
    analysis->setCurMode(curmode);

    if (curmode == AnalysisType::no_An) return MINISPICEOK;

MPI_ENV_0_PID_BEGIN(myid)

    switch (analysis->getCurMode()) {
    case AnalysisType::DC_NR:
        std::cout << "---Curmode is DC_NR" << std::endl;
        break;
    case AnalysisType::DC_Homotopy:
        std::cout << "---Curmode is DC_Homotopy" << std::endl;
        break;
    case AnalysisType::DC_PTran:
        std::cout << "---Curmode is DC_Ptran" << std::endl;
        break;
    case AnalysisType::TRAN:
        std::cout << "---Curmode is Tran" << std::endl;
        break;
    case AnalysisType::no_An:
        std::cout << "---Curmode is no_An" << std::endl;
        break;
    default:
        break;
    }

    std::cout << "---Setup beginning" << std::endl;
    setup(matrix, netlist, analysis, outFileName);
    std::cout << "---Setup ending" << std::endl;

MPI_ENV_0_PID_END

    int error = tranAn(matrix, netlist, analysis, myid);

MPI_ENV_0_PID_BEGIN(myid)
    if (error == MINISPICEERROR) {
        std::cerr << "Tran analysis has faults!" << std::endl;
        local_exit_flag = 1;
        // return MINISPICEERROR;
    }
MPI_ENV_0_PID_END

#ifdef MINISPICESUPERLU
    MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, matrix->superluMatrix->grid.comm);
#else
    MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (global_exit_flag) {
        return MINISPICEERROR;
    }

    unsetup(matrix, netlist, analysis);

    return MINISPICEOK;
    
}

int homotopy(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int myid) {
    return MINISPICEOK;
}

int newton_raphson(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int maxNRIter, int startIdx, int myid) {
MPI_ENV_0_PID_BEGIN(myid)
    // 为 NR 解赋初值
    std::fill(matrix->oldRHS, matrix->oldRHS + matrix->n, 0);
    netlist->loadNodeSet(matrix->oldRHS);

    // std::cout << "------NR solve begin" << std::endl;
MPI_ENV_0_PID_END

    int nr_iters = solve(matrix, netlist, analysis, maxNRIter, startIdx, myid);

MPI_ENV_0_PID_BEGIN(myid)
    // std::cout << "------NR solve end" << std::endl;
MPI_ENV_0_PID_END

    if (nr_iters == MINISPICEERROR || nr_iters > maxNRIter) return NONCONV;

MPI_ENV_0_PID_BEGIN(myid)
    NodeHead &nodeHead = netlist->getNodeHead();
    std::cout << "------";
    for (int i = 0; i < matrix->n; i++) {
        std::cout << "X[" << nodeHead.getID2Name(i + 1) << "]=" << matrix->oldRHS[i] << " ";
    }
    std::cout << std::endl;
MPI_ENV_0_PID_END

    return nr_iters;
}

int pseudoTran(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int myid) {
    double ptran_step_size = analysis->getInitPtranStepSize();
    double max_step = analysis->getMaxPtranStepSize();
    double min_step = analysis->getMinTranStepSize();
    double end_time = analysis->getPtranEndTime();

    int local_exit_flag = 0;
    int global_exit_flag = 0;

    analysis->setPtanStepSize(ptran_step_size);

    int ptran_step_num = 0, ptran_nr_num = 0;
    double ptranTime = 0;
    int ptran_conv = NONCONV;

    // 为 ptran 解赋初值
    std::fill(matrix->oldoldRHS, matrix->oldoldRHS + matrix->n, 0);
    // netlist->loadNodeSet(matrix->oldoldRHS);
    // matrix->oldoldRHS[0] = 4.99922e-05, matrix->oldoldRHS[1] = 0.136364, matrix->oldoldRHS[2] = 2.37682e-18, matrix->oldoldRHS[3] = 9.9984e-05, matrix->oldoldRHS[4] = 4.99922e-11;
    // matrix->oldoldRHS[5] = 9.9985, matrix->oldoldRHS[6] = -0.0149976, matrix->oldoldRHS[7] = -13.6364;
MPI_ENV_0_PID_BEGIN(myid)
    NodeHead &nodeHead = netlist->getNodeHead();
    std::cout << "------initial solution: ";
    for (int i = 0; i < matrix->n; i++) {
        std::cout << "X[" << nodeHead.getID2Name(i + 1) << "]=" << matrix->oldoldRHS[i] << " ";
    }
    std::cout << std::endl;
MPI_ENV_0_PID_END

    while (true) {

        ptranTime += ptran_step_size;

        // 先判断终止结束标志
    MPI_ENV_0_PID_BEGIN(myid)
        if (ptranTime >= end_time || ptran_step_size <= min_step || ptran_step_size >= max_step) {
            std::cerr << "DC ptran failed" << std::endl;
            local_exit_flag = 1;
        }
    MPI_ENV_0_PID_END

    //
    #ifdef MINISPICESUPERLU
        MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, matrix->superluMatrix->grid.comm);
    #else
        MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif

        if (global_exit_flag) {
            return MINISPICEERROR;
        }


        int nr_conv = newton_raphson(matrix, netlist, analysis, analysis->getPtranMaxNRIter(), ptran_step_num, myid);

        ptran_step_num++;
    
    MPI_ENV_0_PID_BEGIN(myid)
        std::cout << "------stepno=" << ptran_step_num << ", timepoint=" << ptranTime << ", pcap=" << analysis->getPCapactor() << ", pind=" << analysis->getPInductor() << ", conv=" << nr_conv << std::endl;
    
        if (nr_conv == NONCONV) {
            // 步长回退1/8
            ptranTime -= ptran_step_size;
            ptran_step_size /= 8;
            analysis->setPtanStepSize(ptran_step_size);
            // std::cout << ptran_step_size << std::endl;
        } else {
            // 步长加2倍
            if (nr_conv < 25)
                ptran_step_size *= 2;
            analysis->setPtanStepSize(ptran_step_size);
            ptran_nr_num += nr_conv;

            double *tmp = matrix->oldoldRHS;
            matrix->oldoldRHS = matrix->oldRHS;
            matrix->oldRHS = tmp;

            // for (int i = 0; i < matrix->n; i++) {
            //     std::cout << "X[" << nodeHead.getID2Name(i + 1) << "]=" << matrix->oldoldRHS[i] << " ";
            // }
            // std::cout << std::endl;

            if (ptConv(matrix, analysis)) {
                ptran_conv = ptran_nr_num;
                local_exit_flag = 1;
                // break;
            }
        }
    MPI_ENV_0_PID_END

        //
        #ifdef MINISPICESUPERLU
            MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, matrix->superluMatrix->grid.comm);
        #else
            MPI_Allreduce(&local_exit_flag, &global_exit_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        #endif

        if (global_exit_flag) {
            break;
        }
        // std::cout << myid << std::endl;
    }

    return ptran_conv;
}

int tranAn(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, int myid) {
    return MINISPICEOK;
}

int trunc();

int ptConv(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis) {
    double rtol = analysis->getPtranRtol(), atol = analysis->getPtranAtol();
    double abs_error = 0.0, rel_error = 0.0;
    for (int i = 0; i < matrix->n; i++) {
        abs_error += (matrix->oldoldRHS[i] - matrix->oldRHS[i]) * (matrix->oldoldRHS[i] - matrix->oldRHS[i]);
    }
    abs_error = std::sqrt(abs_error);
    for (int i = 0; i < matrix->n; i++) {
        rel_error += matrix->oldoldRHS[i] * matrix->oldoldRHS[i];
    }
    rel_error = std::sqrt(rel_error);
    rel_error = abs_error / rel_error;

    if (abs_error <= atol && rel_error <= rtol) return true;
    return false;
}
