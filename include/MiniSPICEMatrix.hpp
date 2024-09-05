
#ifndef MINISPICEMATRIX_H
#define MINISPICEMATRIX_H

// #include <klu.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#ifdef MINISPICEKLU
#include "klu.h"
#endif

#ifdef MINISPICEMUMPS
#include "dmumps_c.h"
#endif

#ifdef MINISPICESUPERLU
#include "superlu_ddefs.h"
#endif

#ifdef MINISPICEEIGEN
#include "Eigen"
#endif

extern const int NA;

struct COOElement {
    int id;
    int row;
    int col;
    int group;
    std::ostringstream JACOutFile;
    COOElement *next;
    bool operator< (const COOElement& b) const {
        return col == b.col ? row < b.row : col < b.col;
    }
    COOElement(int idIn, int rowIn, int colIn):id(idIn), row(rowIn), col(colIn), next(nullptr), group(-1) {}
    COOElement(const COOElement& other):id(other.id), row(other.row), col(other.col), group(other.group), next(nullptr) {
        JACOutFile.str(other.JACOutFile.str());
        JACOutFile.clear();
    }
    COOElement& operator=(const COOElement& other) {
        if (this != &other) {
            id = other.id;
            row = other.row;
            col = other.col;
            group = other.group;
            next = nullptr;
            JACOutFile.str(other.JACOutFile.str());
            JACOutFile.clear(); // 清除可能的错误状态
        }
        return *this;
    }
};

bool compareCOO(const COOElement &a, const COOElement &b);

bool compareCSC(const COOElement &a, const COOElement &b);

bool compareCSR(const COOElement &a, const COOElement &b);

// 在处理雅可比矩阵时添加临时元素
int addCOOElement(COOElement **head, int row, int col, int &elementMapCount);


#ifdef MINISPICEKLU
struct KLUMatrix {
    klu_common KLUmatrixCommon;
    klu_symbolic *KLUmatrixSymbolic;
    klu_numeric *KLUmatrixNumeric;
    int n;
    int *KLUmatrixAp;
    int *KLUmatrixAi;
    double *KLUmatrixAx;
    KLUMatrix():n(0), KLUmatrixSymbolic(nullptr), KLUmatrixNumeric(nullptr), KLUmatrixAp(nullptr), KLUmatrixAi(nullptr), KLUmatrixAx(nullptr) {}
    // ~KLUMatrix() {
    //     klu_free_symbolic(&KLUmatrixSymbolic, &KLUmatrixCommon);
    //     klu_free_numeric(&KLUmatrixNumeric, &KLUmatrixCommon);
    //     delete KLUmatrixAp;
    //     delete KLUmatrixAi;
    //     delete KLUmatrixAx;
    //     KLUmatrixAp = KLUmatrixAi = nullptr;
    //     KLUmatrixAx = nullptr;
    //     n = 0;
    // }
    void destroy() {
        klu_free_symbolic(&KLUmatrixSymbolic, &KLUmatrixCommon);
        klu_free_numeric(&KLUmatrixNumeric, &KLUmatrixCommon);
        delete KLUmatrixAp;
        delete KLUmatrixAi;
        delete KLUmatrixAx;
        KLUmatrixAp = KLUmatrixAi = nullptr;
        KLUmatrixAx = nullptr;
        n = 0;
    }
};
#endif

#ifdef MINISPICEMUMPS
struct MUMPSMatrix {
    #define JOB_INIT -1
    #define JOB_END -2
    #define USE_COMM_WORLD -987654
    #define ICNTL(I) icntl[(I)-1] //与mumps的文档对齐
    DMUMPS_STRUC_C id;
    int myid;
    int n;
    int nnz;
    int *rowidx;
    int *colidx;
    double *elementVal;
    MUMPSMatrix():rowidx(nullptr), colidx(nullptr), elementVal(nullptr), n(0), nnz(0), myid(0) {}
    // ~MUMPSMatrix() {
    //     id.job = JOB_END;
    //     dmumps_c(&id);
    //     if (myid == 0) {
    //         delete rowidx;
    //         delete colidx;
    //         delete elementVal;
    //         rowidx = colidx = nullptr;
    //         elementVal = nullptr;
    //     }
    //     n = nnz = 0;
    //     // myid 将一直生效
    // }
    void destroy() {
        id.job = JOB_END;
        dmumps_c(&id);
        if (myid == 0) {
            delete rowidx;
            delete colidx;
            delete elementVal;
            rowidx = colidx = nullptr;
            elementVal = nullptr;
        }
        n = nnz = 0;
        // myid 将一直生效
    }
};
#endif

#ifdef MINISPICESUPERLU
struct SuperLUMatrix {
    superlu_dist_options_t options;  // solver参数控制选项
    SuperLUStat_t stat;
    SuperMatrix A;
    dScalePermstruct_t ScalePermstruct;
    dLUstruct_t LUstruct;
    dSOLVEstruct_t SOLVEstruct;
    gridinfo_t grid;
    int      omp_mpi_level;
    int      world_size;
    double   berr[1];
    
    // 局部参与运算的矩阵[所有进程各自的行快]
    int      local_m, local_n;
    int      local_nnz;
    int      fst_row;
    double   *local_CSRValue;
    int      *local_rowptr, *local_rowptr_first;
    int      *local_colind, *local_colind_first;
    double   *local_rhs;

    // 全局的矩阵的值，生命周期随矩阵生命周期【只有0进程有所有信息】
    int      global_m, global_n;
    int      global_nnz;
    int      *global_rowptr;
    int      *global_colind;
    double   *global_CSRValue;
    /*! 它的原型在上层结构中 !*/
    double   *global_rhs;

    // 分发相关，仅在进程0上定义发送计数和偏移量
    int *sendcounts_csr_colind_val, *sendcounts_csr_rowptr, *sendcounts_csr_rhs;
    int *displs_csr_colind_val, *displs_csr_rowptr, *displs_csr_rhs;

    SuperLUMatrix() {
        sendcounts_csr_colind_val = sendcounts_csr_rowptr = sendcounts_csr_rhs = nullptr;
        displs_csr_colind_val = displs_csr_rowptr = displs_csr_rhs = nullptr;
        global_colind = global_rowptr = nullptr;
        local_colind = local_rowptr = local_rowptr_first = local_colind_first = nullptr;
        global_rhs = local_rhs = local_CSRValue = global_CSRValue = nullptr;
        // berr = nullptr;

        global_m = global_n = global_nnz = local_n = local_m = local_nnz = fst_row = -1;
    }

    // ~SuperLUMatrix() {
    //     int iam = grid.iam;
    //     if (iam == 0) {
    //         delete global_rowptr;
    //         delete global_colind;
    //         delete global_CSRValue;
            
    //         delete sendcounts_csr_colind_val;
    //         delete sendcounts_csr_rowptr;
    //         delete sendcounts_csr_rhs;
    //         delete displs_csr_colind_val;
    //         delete displs_csr_rowptr;
    //         delete displs_csr_rhs;
    //     }
    //     // delete local_CSRValue;
    //     // delete local_rowptr;
    //     // delete local_colind;
    //     delete local_rhs;

    //     sendcounts_csr_colind_val = sendcounts_csr_rowptr = sendcounts_csr_rhs = nullptr;
    //     displs_csr_colind_val = displs_csr_rowptr = displs_csr_rhs = nullptr;
    //     global_colind = global_rowptr = nullptr;
    //     // local_rhs = local_colind = local_rowptr = local_CSRValue = nullptr;
    //     global_CSRValue = global_rhs = local_rhs = nullptr;

    //     global_m = global_n = global_nnz = local_n = local_m = local_nnz = fst_row = -1;


    //     // PStatFree(&stat);
    //     Destroy_CompRowLoc_Matrix_dist(&A);
    //     // dDestroy_LU(local_n, &grid, &LUstruct);
    //     dScalePermstructFree(&ScalePermstruct);
    //     dLUstructFree(&LUstruct);         /* Deallocate the structure of L and U.*/
    //     if ( options.SolveInitialized ) {
    //         dSolveFinalize(&options, &SOLVEstruct);
    //     }
    //     SUPERLU_FREE(berr);
    // }
    void destroy() {
        int iam = grid.iam;
        if (iam == 0) {
            delete global_rowptr;
            delete global_colind;
            delete global_CSRValue;

            delete local_colind_first;
            delete local_rowptr_first;
            
            delete sendcounts_csr_colind_val;
            delete sendcounts_csr_rowptr;
            delete sendcounts_csr_rhs;
            delete displs_csr_colind_val;
            delete displs_csr_rowptr;
            delete displs_csr_rhs;
        }
        // delete local_CSRValue;
        // delete local_rowptr;
        // delete local_colind;
        delete local_rhs;

        sendcounts_csr_colind_val = sendcounts_csr_rowptr = sendcounts_csr_rhs = nullptr;
        displs_csr_colind_val = displs_csr_rowptr = displs_csr_rhs = nullptr;
        global_colind = global_rowptr = nullptr;
        // local_rhs = local_colind = local_rowptr = local_CSRValue = nullptr;
        local_rowptr_first = local_colind_first = nullptr;
        global_CSRValue = global_rhs = local_rhs = nullptr;

        global_m = global_n = global_nnz = local_n = local_m = local_nnz = fst_row = -1;


        PStatFree(&stat);
        Destroy_CompRowLoc_Matrix_dist(&A);
        // dDestroy_LU(local_n, &grid, &LUstruct);
        dScalePermstructFree(&ScalePermstruct);
        dLUstructFree(&LUstruct);         /* Deallocate the structure of L and U.*/
        if ( options.SolveInitialized ) {
            dSolveFinalize(&options, &SOLVEstruct);
        }
        // SUPERLU_FREE(berr);
    }

};

// 矩阵只在外层创建一次
// 传入的my_martix必须已经构建好了全局的CSR矩阵
int dcreate_matrix(SuperLUMatrix *my_matrix, int startIdx);

#endif

#ifdef MINISPICEEIGEN
struct EigenMatrix {
    
};
#endif

// TODO
struct MiniSPICEMatrix {
    void destroy();
};


struct SMPMatrix {
#if defined(MINISPICEKLU)
    KLUMatrix *kluMatrix;
#elif defined(MINISPICEMUMPS)
    MUMPSMatrix *mumpsMatrix;
#elif defined(MINISPICESUPERLU)
    SuperLUMatrix *superluMatrix;
#elif defined(MINISPICEEIGEN)
    EigenMatrix *eigenMatrix;
#else
    MiniSPICEMatrix *miniSPICEMatrix;
#endif
    
    //
    double *RHS;
    double *oldRHS;
    double *oldoldRHS; // DC 下， oldRHS == oldoldRHS
    double **elementMap;
    int elementNNZ;
    int elementMapCount;
    int elementMapCountLocal;
    int n;

    //
    COOElement *coo;

    //
    SMPMatrix():RHS(nullptr), oldRHS(nullptr), oldoldRHS(nullptr), elementMap(nullptr), elementNNZ(0), elementMapCount(0), elementMapCountLocal(0), n(0), coo(nullptr){
    #if defined(MINISPICEKLU)
        kluMatrix = new KLUMatrix();
    #elif defined(MINISPICEMUMPS)
        mumpsMatrix = new MUMPSMatrix();
    #elif defined(MINISPICESUPERLU)
        superluMatrix = new SuperLUMatrix();
    #elif defined(MINISPICEEIGEN)
        eigenMatrix = new EigenMatrix();
    #else
        miniSPICEMatrix = new MiniSPICEMatrix;
    #endif
    }

    void destroySMPMatrix() {
        if (RHS != nullptr) delete RHS;
        if (oldRHS != nullptr) delete oldRHS;
        if (oldoldRHS != nullptr) delete oldoldRHS;
        if (elementMap != nullptr) delete elementMap;
        RHS = oldoldRHS = oldRHS = nullptr;
        elementMap = nullptr;
        coo = nullptr;
        n = elementNNZ = elementMapCount = elementMapCountLocal = 0;

    #if defined(MINISPICEKLU)
        // delete kluMatrix;
        kluMatrix->destroy();
    #elif defined(MINISPICEMUMPS)
        // delete mumpsMatrix;
        mumpsMatrix->destroy();
    #elif defined(MINISPICESUPERLU)
        // delete superluMatrix;
        superluMatrix->destroy();
    #elif defined(MINISPICEEIGEN)
        // delete eigenMatrix;
        eigenMatrix->destroy();
    #else
        // delete miniSPICEMatrix;
        miniSPICEMatrix->destroy();
    #endif
    }

    ~SMPMatrix() {
        #if defined(MINISPICEKLU)
            delete kluMatrix;
            // kluMatrix->destroy();
        #elif defined(MINISPICEMUMPS)
            delete mumpsMatrix;
            // mumpsMatrix->destroy();
        #elif defined(MINISPICESUPERLU)
            delete superluMatrix;
            // superluMatrix->destroy();
        #elif defined(MINISPICEEIGEN)
            delete eigenMatrix;
            // eigenMatrix->destroy();
        #else
            delete miniSPICEMatrix;
            // miniSPICEMatrix->destroy();
        #endif
    }

};

#endif
