#include<algorithm>
#include "minispice.hpp"

int setup(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, std::string outFileName) {

    AnalysisType curMode = analysis->getCurMode();
    if (curMode == AnalysisType::no_An) return MINISPICEOK;
    if (curMode == AnalysisType::TRAN ) outFileName = outFileName + "_tran";
    else if (curMode == AnalysisType::DC_PTran) outFileName = outFileName + "_ptran";
    else outFileName = outFileName + "_op";

    std::cout << "------Outfile is " << outFileName << std::endl;

    std::vector<std::ostringstream> &FoutFile = analysis->getFOutFileAndClear();

    std::cout << "------Device list scanning begin" << std::endl;
    
    // 1. setup主要阶段，对FX、jac矩阵进行信息处理
    NodeHead &nodeHead = netlist->getNodeHead();
    Node *nodePtr = nodeHead.getNode();
    while (nodePtr != nullptr) {
        if (nodePtr->getNodeNum() != 0) {
            Connections *conList = nodePtr->getConList();
            bool bjtflag = false;
            while (conList != nullptr) {
                Device *device = conList->device;
                if (device->getType() == DeviceType::BJT) bjtflag = true;
                // std::cout << device->getName() << " " << matrix->elementMapCount << " " << device->getNodeNum(conList->conNum) << " " << conList->conNum << std::endl;
                int error = device->setup(&(matrix->coo), matrix->elementMapCount, FoutFile, analysis, conList->conNum);
                if (error == MINISPICEERROR) {
                    std::cerr << "setup phase has error occurred!" << std::endl;
                    return MINISPICEERROR;
                }
                conList = conList->next;
            }
            // 处理ptran
            if (curMode == AnalysisType::DC_PTran && bjtflag) {
                int t = nodePtr->getNodeNum() - 1;
                // FX
                // matrix->RHS[t] += analysis->getPCapactor() / analysis->getPTranStepSize() * (matrix->oldRHS[t] - matrix->oldoldRHS[t]);
                FoutFile[t] << " + C / h * (X_" << nodePtr->getNodeName() << " - X'_" << nodePtr->getNodeName() << ")";
                //JAC
                // *matrix->elementMap[matrix->elementMapCount++] += analysis->getPCapactor() / analysis->getPTranStepSize();
                if (addCOOElement(&(matrix->coo), t, t, matrix->elementMapCount) == NA) return MINISPICEERROR;
                matrix->coo->JACOutFile << " + C / h";
            }
        }
        nodePtr = nodePtr->getNext();
    }

    std::cout << "------Device list scanning end" << std::endl;
    std::cout << "------Matrix generation start" << std::endl;
    std::cout << "---------Matrix mapping process start" << std::endl;


    // 2. 处理矩阵存储
    matrix->n = FoutFile.size();
    matrix->RHS = new double[matrix->n];
    matrix->oldRHS = new double[matrix->n];
    matrix->oldoldRHS = new double[matrix->n];
    matrix->elementMap = new double*[matrix->elementMapCount];

    // 2.1 处理映射表
    std::vector<COOElement> coos_vec;
    COOElement *coo_ptr = matrix->coo, *temp;
    while (coo_ptr != nullptr) {
        temp = coo_ptr->next;
        coos_vec.push_back(*coo_ptr);
        delete coo_ptr;
        coo_ptr = temp;
    }
    matrix->coo = nullptr;

    if (matrix->elementMapCount != coos_vec.size()) {
        std::cerr << "coo -> coo_vec failed!" << std::endl;
        return MINISPICEERROR;
    }

    // 2.2 按不同类型矩阵存储处理
#if defined(MINISPICEKLU)
    std::sort(coos_vec.begin(), coos_vec.end(), compareCSC);
#elif defined(MINISPICESUPERLU)
    std::sort(coos_vec.begin(), coos_vec.end(), compareCSR);
#else
    std::sort(coos_vec.begin(), coos_vec.end(), compareCOO);
#endif

    // 2.3 仅从元素 “去重”
    if (coos_vec.size() != 0) coos_vec[0].group = matrix->elementNNZ++;
    // std::cout << coos_vec[0].group << " " << coos_vec[0].id  << " " << coos_vec[0].row << " " << coos_vec[0].col << std::endl;
    for (int i = 1; i < coos_vec.size(); i++) {
        if (coos_vec[i].row == coos_vec[i - 1].row && coos_vec[i].col == coos_vec[i - 1].col)
            coos_vec[i].group = coos_vec[i - 1].group;
        else 
            coos_vec[i].group = matrix->elementNNZ++;
        // std::cout << coos_vec[i].group << " " << coos_vec[i].id  << " " << coos_vec[i].row << " " << coos_vec[i].col << std::endl;
    }

    std::cout << "---------Matrix mapping process end" << std::endl;

    // 2.4 矩阵空间申请及对应map
#if defined(MINISPICEKLU)
    std::cout << "---------Matrix klu process start" << std::endl;
    matrix->kluMatrix->n = matrix->n;
    matrix->kluMatrix->KLUmatrixAp = new int[matrix->n + 1];
    matrix->kluMatrix->KLUmatrixAp[0] = 0;
    matrix->kluMatrix->KLUmatrixAx = new double[matrix->elementNNZ];
    matrix->kluMatrix->KLUmatrixAi = new int[matrix->elementNNZ];

    for (int i = 0; i < coos_vec.size(); i++) {
         // 记录非零元的映射关系
        matrix->elementMap[coos_vec[i].id] = &matrix->kluMatrix->KLUmatrixAx[coos_vec[i].group];
        // std::cout << matrix->elementMap[coos_vec[i].id] << " ";
        // 记录非零元的行号
        matrix->kluMatrix->KLUmatrixAi[coos_vec[i].group] = coos_vec[i].row;
        // 记录前一列的非零元数目
        if (i != 0 && coos_vec[i].col != coos_vec[i - 1].col) {
            matrix->kluMatrix->KLUmatrixAp[coos_vec[i].col] = coos_vec[i].group;
            // std::cout << coos_vec[i].group << std::endl;
        }
    }
    // std::cout << std::endl;
    matrix->kluMatrix->KLUmatrixAp[matrix->n] = matrix->elementNNZ;
    // std::cout << matrix->elementNNZ << std::endl;
    // for (int i = 0; i < matrix->n + 1; i++) {
    //     std::cout << matrix->kluMatrix->KLUmatrixAp[i] << " ";
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < matrix->elementNNZ; i++) {
    //     std::cout << matrix->kluMatrix->KLUmatrixAi[i] << " ";
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < matrix->elementNNZ; i++) {
    //     std::cout << matrix->kluMatrix->KLUmatrixAx[i] << " ";
    // }
    // std::cout << std::endl;

    std::cout << "---------Matrix klu process end" << std::endl;

#elif defined(MINISPICESUPERLU)
    std::cout << "---------Matrix superlu process start" << std::endl;
    matrix->superluMatrix->global_m = matrix->n;
    matrix->superluMatrix->global_n = matrix->n;
    matrix->superluMatrix->global_nnz = matrix->elementNNZ;
    matrix->superluMatrix->global_rowptr = new int[matrix->superluMatrix->global_m + 1];
    matrix->superluMatrix->global_colind = new int[matrix->superluMatrix->global_nnz];
    matrix->superluMatrix->global_CSRValue = new double[matrix->superluMatrix->global_nnz];
    matrix->superluMatrix->global_rhs = matrix->RHS;

    for (int i = 0; i < coos_vec.size(); i++) {
         // 记录非零元的映射关系
        matrix->elementMap[coos_vec[i].id] = &matrix->superluMatrix->global_CSRValue[coos_vec[i].group];
        // 记录非零元的列号
        matrix->superluMatrix->global_colind[coos_vec[i].group] = coos_vec[i].col;
        // 记录前一列的非零元数目
        if (i != 0 && coos_vec[i].row != coos_vec[i - 1].row) {
            matrix->superluMatrix->global_rowptr[coos_vec[i].row] = coos_vec[i].group;
        }
    }
    matrix->superluMatrix->global_rowptr[0] = 0;
    matrix->superluMatrix->global_rowptr[matrix->superluMatrix->global_m] = matrix->elementNNZ;

    // std::cout << matrix->elementNNZ << std::endl;
    // for (int i = 0; i < matrix->n + 1; i++) {
    //     std::cout << matrix->superluMatrix->global_rowptr[i] << " ";
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < matrix->elementNNZ; i++) {
    //     std::cout << matrix->superluMatrix->global_colind[i] << " ";
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < matrix->elementNNZ; i++) {
    //     std::cout << matrix->superluMatrix->global_CSRValue[i] << " ";
    // }
    // std::cout << std::endl;

    std::cout << "---------Matrix superlu process end" << std::endl;

#elif defined(MINISPICEMUMPS)
    std::cout << "---------Matrix mumps process start" << std::endl;
    matrix->mumpsMatrix->n = matrix->n;
    matrix->mumpsMatrix->nnz = matrix->elementNNZ;
    matrix->superluMatrix->rowidx = new int[matrix->mumpsMatrix->nnz];
    matrix->superluMatrix->colidx = new int[matrix->mumpsMatrix->nnz];
    matrix->superluMatrix->elementVal = new double[matrix->mumpsMatrix->nnz];

    for (int i = 0; i < coos_vec.size(); i++) {
         // 记录非零元的映射关系
        matrix->elementMap[coos_vec[i].id] = &matrix->superluMatrix->elementVal[coos_vec[i].group];
        // 记录非零元的行列号
        matrix->superluMatrix->rowidx[coos_vec[i].group] = coos_vec[i].row;
        matrix->superluMatrix->colidx[coos_vec[i].group] = coos_vec[i].col;
    }
    std::cout << "---------Matrix mumps process end" << std::endl;

#endif

    std::cout << "------Matrix generation end" << std::endl;
    std::cout << "------Circuit equation(kcl/jac) generation start--" << outFileName + "_FX_.txt" << std::endl;
    
    // 4. 输出信息到文件
    std::ofstream outfile1(outFileName + "_FX_.txt");
    if (outfile1.is_open()) {
        for (int i = 0; i < FoutFile.size(); i++) {
            outfile1 << "F_" << nodeHead.getID2Name(i + 1) << " =" << FoutFile[i].str() << std::endl;
        }
        outfile1.close();
    } else {
        std::cerr << "Failed to open the file." << std::endl;
        return MINISPICEERROR;
    }
    std::ofstream outfile2(outFileName + "_JAC_.txt");
    if (outfile2.is_open()) {
        for (int i = 0; i < coos_vec.size(); i++) {
            if (i != 0 && coos_vec[i].row == coos_vec[i - 1].row && coos_vec[i].col == coos_vec[i - 1].col)
                outfile2 << coos_vec[i].JACOutFile.str();
            else {
                if (i != 0) outfile2 << std::endl;
                outfile2 << "JAC(" << nodeHead.getID2Name(coos_vec[i].row + 1) << "," << nodeHead.getID2Name(coos_vec[i].col + 1)
                         << ") =" << coos_vec[i].JACOutFile.str();
            }
        }
        outfile2.close();
    } else {
        std::cerr << "Failed to open the file." << std::endl;
        return MINISPICEERROR;
    }
    std::cout << "------Circuit equation(kcl/jac) generation end--" << outFileName + "_JAC_.txt" << std::endl;

    return MINISPICEOK;
}


int unsetup(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis) {

    matrix->destroySMPMatrix();

    return MINISPICEOK;
}


int load(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis) {

    AnalysisType curMode = analysis->getCurMode();
    NodeHead &nodeHead = netlist->getNodeHead();
    Node *nodePtr = nodeHead.getNode();

    matrix->elementMapCountLocal = 0;
    for (int i = 0; i < matrix->elementMapCount; i++) {
        *matrix->elementMap[i] = 0;
    }
    // for (int i = 0; i < matrix->elementNNZ; i++) {
    //     std::cout << matrix->kluMatrix->KLUmatrixAx[i] << " ";
    // }
    // std::cout << std::endl;
    for (int i = 0; i < matrix->n; i++) {
        matrix->RHS[i] = 0;
    }
    // netlist->loadNodeSet(matrix->oldRHS);

    while (nodePtr != nullptr) {
        if (nodePtr->getNodeNum() != 0) {
            Connections *conList = nodePtr->getConList();
            bool bjtflag = false;
            while (conList != nullptr) {
                Device *device = conList->device;
                if (device->getType() == DeviceType::BJT) bjtflag = true;
                // std::cout << device->getName() << " " << matrix->elementMapCountLocal << " " << device->getNodeNum(conList->conNum) << " " << conList->conNum << std::endl;
                int error = device->genKCLEquation(matrix, analysis, conList->conNum);
                if (error == MINISPICEERROR) {
                    std::cerr << "load phase has error occurred!" << std::endl;
                    return MINISPICEERROR;
                }
                error = device->genKCLJAC(matrix, analysis, conList->conNum);
                if (error == MINISPICEERROR) {
                    std::cerr << "load phase has error occurred!" << std::endl;
                    return MINISPICEERROR;
                }
                conList = conList->next;
            }
            // 处理ptran
            if (curMode == AnalysisType::DC_PTran && bjtflag) {
                int t = nodePtr->getNodeNum() - 1;
                // FX
                matrix->RHS[t] += analysis->getPCapactor() / analysis->getPTranStepSize() * (matrix->oldRHS[t] - matrix->oldoldRHS[t]);
                // FoutFile[t] << " + C / h * (X_" << nodePtr->getNodeName() << " - X'_" << nodePtr->getNodeName() << ")";
                //JAC
                *matrix->elementMap[matrix->elementMapCountLocal++] += analysis->getPCapactor() / analysis->getPTranStepSize();
                // if (addCOOElement(&(matrix->coo), t, t, matrix->elementMapCount) == NA) return MINISPICEERROR;
                // matrix->coo->JACOutFile << " + C / h";
            }
        }
        nodePtr = nodePtr->getNext();
    }
    if (matrix->elementMapCountLocal != matrix->elementMapCount) {
        std::cerr << "load phase has error occurred [matrix->elementMapCountLocal != matrix->elementMapCount]!" << std::endl;
        std::cerr << matrix->elementMapCountLocal << " " << matrix->elementMapCount << std::endl;
        return MINISPICEERROR;
    }
    // std::cout << matrix->elementMapCountLocal << " " << matrix->elementMapCount << " " << analysis->getPTranStepSize() << " " << analysis->getPCapactor() << " " << analysis->getPInductor() << std::endl;
    // for (int i = 0; i < matrix->elementNNZ; i++) {
    //     std::cout << matrix->kluMatrix->KLUmatrixAx[i] << " ";
    //     // std::cout << matrix->superluMatrix->global_CSRValue[i] << " ";
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < matrix->n; i++) {
    //     std::cout << matrix->RHS[i] << " ";
    // }
    // std::cout << std::endl;
    // return MINISPICEERROR;
    return MINISPICEOK;
}
