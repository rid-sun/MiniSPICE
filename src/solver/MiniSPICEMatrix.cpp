#include "minispice.hpp"
#include "MiniSPICEMatrix.hpp"

bool compareCOO(const COOElement &a, const COOElement &b) {
    if (a.row == b.row && a.col == b.col) {
        return a.id < b.id;
    } else if (a.row == b.row) {
        return a.col < b.col;
    } else {
        return a.row < b.row;
    }
}

bool compareCSC(const COOElement &a, const COOElement &b) {
    if (a.row == b.row && a.col == b.col) {
        return a.id < b.id;
    } else if (a.col == b.col) {
        return a.row < b.row;
    } else {
        return a.col < b.col;
    }
}

bool compareCSR(const COOElement &a, const COOElement &b) {
    if (a.row == b.row && a.col == b.col) {
        return a.id < b.id;
    } else if (a.row == b.row) {
        return a.col < b.col;
    } else {
        return a.row < b.row;
    }
}

// 在处理雅可比矩阵时添加临时元素
int addCOOElement(COOElement **head, int row, int col, int &elementMapCount) {
    int id = elementMapCount++;
    COOElement *ptr = new COOElement(id, row, col);
    if (ptr == nullptr) return NA;
    ptr->next = *head;
    *head = ptr;
    return 0;
}



#ifdef MINISPICESUPERLU

// 矩阵只在外层创建一次
// 传入的my_martix必须已经构建好了全局的CSR矩阵
int dcreate_matrix(SuperLUMatrix *my_matrix, int startIdx) {
    int_t    m, n;
    int_t    m_loc, m_loc_last;
    int      iam;

    SuperMatrix *A = &my_matrix->A;
    gridinfo_t *grid = &my_matrix->grid;

    iam = grid->iam;

    int     *msg_row;             /* 分别为行数、非零元数、行偏移量*/
    int     *msg_nnz;             /* 分别为行数、非零元数、行偏移量*/
    int     *msg_offset_row;      /* 分别为行数、非零元数、行偏移量*/

if (startIdx == 1){

    if (!iam) {
        
        m = my_matrix->global_m;
        n = my_matrix->global_n;

        m_loc = m / (grid->nprow * grid->npcol);
        int offset;
        int num_process = grid->nprow * grid->npcol;

        msg_row = new int[num_process];
        msg_nnz = new int[num_process];
        msg_offset_row = new int[num_process];


        // 计算分发块的大小-RHS
        my_matrix->displs_csr_rhs = new int[num_process];
        my_matrix->sendcounts_csr_rhs = new int[num_process];
        offset = 0;
        for (int i = 0; i < num_process - 1; i++) {
            msg_row[i] = m_loc;
            my_matrix->sendcounts_csr_rhs[i] = m_loc;
            my_matrix->displs_csr_rhs[i] = offset;
            msg_offset_row[i] = offset;
            offset += my_matrix->sendcounts_csr_rhs[i];
        }
        my_matrix->sendcounts_csr_rhs[num_process - 1] = my_matrix->global_m - m_loc * (num_process - 1);
        my_matrix->displs_csr_rhs[num_process - 1] = offset;
        msg_offset_row[num_process - 1] = offset;
        msg_row[num_process - 1] = my_matrix->global_m - m_loc * (num_process - 1);
        // if ((m_loc * num_process) != m) {
        //     m_loc_last = m - m_loc * (num_process - 1);
        //     my_matrix->sendcounts_csr_rhs[num_process - 1] = m_loc_last;
        // }

        // // 计算分发块的大小-colind
        // my_matrix->displs_csr_colind = new int[num_process];
        // my_matrix->sendcounts_csr_colind = new int[num_process];
        // offset = 0;
        // for (int i = 0; i < num_process - 1; i++) {
        //     my_matrix->sendcounts_csr_colind[i] = my_matrix->global_colind[m_loc * (i + 1)] - my_matrix->global_colind[m_loc * i];
        //     my_matrix->displs_csr_colind[i] = offset;
        //     offset += my_matrix->sendcounts_csr_colind[i];
        // }
        // my_matrix->sendcounts_csr_colind[num_process - 1] = my_matrix->global_nnz - offset;
        // my_matrix->displs_csr_colind[i] = offset;
        // // if ((m_loc * num_process) != m) {
        // //     m_loc_last = m - m_loc * (num_process - 1);
        // //     my_matrix->sendcounts_csr_colind[num_process - 1] = m_loc_last;
        // // }

        // 计算分发块的大小-rowptr
        my_matrix->displs_csr_rowptr = new int[num_process];
        my_matrix->sendcounts_csr_rowptr = new int[num_process];
        offset = 0;
        for (int i = 0; i < num_process - 1; i++) {
            my_matrix->sendcounts_csr_rowptr[i] = m_loc + 1;
            my_matrix->displs_csr_rowptr[i] = offset;
            offset += m_loc;
        }
        my_matrix->sendcounts_csr_rowptr[num_process - 1] = my_matrix->global_m - m_loc * (num_process - 1) + 1;
        my_matrix->displs_csr_rowptr[num_process - 1] = offset;

        // // 计算分发块的大小-csrval
        // my_matrix->displs_csr_val = new int[num_process];
        // my_matrix->sendcounts_csr_val = new int[num_process];
        // offset = 0;
        // for (int i = 0; i < num_process - 1; i++) {
        //     msg_nnz[i] = my_matrix->sendcounts_csr_val[i] = my_matrix->global_CSRValue[m_loc * (i + 1)] - my_matrix->global_CSRValue[m_loc * i];
        //     my_matrix->displs_csr_val[i] = offset;
        //     offset += my_matrix->sendcounts_csr_val[i];
        // }
        // msg_nnz[num_process - 1] = my_matrix->sendcounts_csr_val[num_process - 1] = my_matrix->global_nnz - offset;
        // my_matrix->displs_csr_val[i] = offset;

        // 计算分发块的大小-colind_csrval
        my_matrix->displs_csr_colind_val = new int[num_process];
        my_matrix->sendcounts_csr_colind_val = new int[num_process];
        offset = 0;
        for (int i = 0; i < num_process - 1; i++) {
            msg_nnz[i] = my_matrix->sendcounts_csr_colind_val[i] = my_matrix->global_rowptr[m_loc * (i + 1)] - my_matrix->global_rowptr[m_loc * i];
            my_matrix->displs_csr_colind_val[i] = offset;
            offset += my_matrix->sendcounts_csr_colind_val[i];
        }
        msg_nnz[num_process - 1] = my_matrix->sendcounts_csr_colind_val[num_process - 1] = my_matrix->global_nnz - offset;
        my_matrix->displs_csr_colind_val[num_process - 1] = offset;

    }

     // 广播 m n
    MPI_Bcast(&m, 1, MPI_INT, 0, grid->comm);
    MPI_Bcast(&n, 1, MPI_INT, 0, grid->comm);

    // 分发 local_m
    MPI_Scatter(
        msg_row,                      // 发送缓冲区，仅在root进程有意义
        1,                            // 每个进程要接收的数据数量
        MPI_INT,                      // 数据类型
        &my_matrix->local_m,          // 接收缓冲区，每个进程都会有
        1,                            // 每个进程要接收的数据数量
        MPI_INT,                      // 数据类型
        0,                            // root进程的rank
        grid->comm                    // 通信器
    );

    // 分发 local_nnz
    MPI_Scatter(msg_nnz, 1, MPI_INT, &my_matrix->local_nnz, 1, MPI_INT, 0, grid->comm);

    // 分发 fst_row
    MPI_Scatter(msg_offset_row, 1, MPI_INT, &my_matrix->fst_row, 1, MPI_INT, 0, grid->comm);
    
    // 创建接受缓冲区
    my_matrix->local_CSRValue = new double[my_matrix->local_nnz];
    my_matrix->local_rowptr = new int[my_matrix->local_m + 1];
    my_matrix->local_rowptr_first = new int[my_matrix->local_m + 1];
    my_matrix->local_colind = new int[my_matrix->local_nnz];
    my_matrix->local_colind_first = new int[my_matrix->local_nnz];
    my_matrix->local_rhs = new double[my_matrix->local_m];

    // 分发 local_CSRValue
    MPI_Scatterv(
        my_matrix->global_CSRValue,           // 发送缓冲区（仅对root进程有意义）
        my_matrix->sendcounts_csr_colind_val, // 各进程接收的数据大小（仅在root进程有意义）
        my_matrix->displs_csr_colind_val,     // 各进程接收数据的位移（仅在root进程有意义）
        MPI_DOUBLE,                           // 数据类型
        my_matrix->local_CSRValue,            // 接收缓冲区
        my_matrix->local_nnz,                 // 当前进程接收的数据大小
        MPI_DOUBLE,                           // 数据类型
        0,                                    // root进程的rank
        grid->comm                            // 通信器
    );

    // 分发 local_colind
    MPI_Scatterv(my_matrix->global_colind, my_matrix->sendcounts_csr_colind_val, my_matrix->displs_csr_colind_val, 
                MPI_INT, my_matrix->local_colind, my_matrix->local_nnz, MPI_INT, 0, grid->comm);

    // 分发 local_rowptr
    MPI_Scatterv(my_matrix->global_rowptr, my_matrix->sendcounts_csr_rowptr, my_matrix->displs_csr_rowptr, 
                MPI_INT, my_matrix->local_rowptr, my_matrix->local_m + 1, MPI_INT, 0, grid->comm);
    
    // 分发 local_rhs
    MPI_Scatterv(my_matrix->global_rhs, my_matrix->sendcounts_csr_rhs, my_matrix->displs_csr_rhs, 
                MPI_DOUBLE, my_matrix->local_rhs, my_matrix->local_m, MPI_DOUBLE, 0, grid->comm);
    
    // 后处理，处理rowptr坐标
    for (int i = 1; i < my_matrix->local_m + 1; i++) {
        my_matrix->local_rowptr[i] -= my_matrix->local_rowptr[0];
    }
    my_matrix->local_rowptr[0] = 0;

    /* 记录备份local_colind和local_rowptr，因为在执行LU等计算之后，它的排列顺序发生改变了 */
    memcpy(my_matrix->local_colind_first, my_matrix->local_colind, my_matrix->local_nnz * sizeof(int));
    memcpy(my_matrix->local_rowptr_first, my_matrix->local_rowptr, (my_matrix->local_m + 1) * sizeof(int));

    // 创建矩阵
    /* Set up the local A in NR_loc format */
    dCreate_CompRowLoc_Matrix_dist(A, m, n, my_matrix->local_nnz, my_matrix->local_m, my_matrix->fst_row,
				   my_matrix->local_CSRValue, my_matrix->local_colind, my_matrix->local_rowptr,
				   SLU_NR_loc, SLU_D, SLU_GE);
    

    // 释放空间
    if (iam == 0) {
        delete msg_row;
        delete msg_nnz;
        delete msg_offset_row;
    }
    
} else {

    // 只需要复制矩阵元素和有端项即可
    // 分发 local_rhs
    MPI_Scatterv(my_matrix->global_rhs, my_matrix->sendcounts_csr_rhs, my_matrix->displs_csr_rhs, 
                MPI_DOUBLE, my_matrix->local_rhs, my_matrix->local_m, MPI_DOUBLE, 0, grid->comm);
    // 分发 local_CSRValue
    MPI_Scatterv(my_matrix->global_CSRValue, my_matrix->sendcounts_csr_colind_val, my_matrix->displs_csr_colind_val, 
                MPI_DOUBLE, my_matrix->local_CSRValue, my_matrix->local_nnz, MPI_DOUBLE, 0, grid->comm);

    // 由于元素排列发生改变，需要重新赋值
    memcpy(my_matrix->local_colind, my_matrix->local_colind_first, my_matrix->local_nnz * sizeof(int));
    memcpy(my_matrix->local_rowptr, my_matrix->local_rowptr_first, (my_matrix->local_m + 1) * sizeof(int));

}
    return MINISPICEOK;
}

#endif
