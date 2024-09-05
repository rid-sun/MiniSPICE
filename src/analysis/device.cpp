#include "device.hpp"

// 常量声明
// const int NameLength = 80;
const int BufLength = 3000;
const int NA = -1;
const int OK = 0;


/****************************************************************************************************/
/***********************************   Device 类定义           ***************************************/
/****************************************************************************************************/
/****************************************************************************************************/

void Device::connect(int conNum, Node *nodeIn) {
    assert(conNum >= 0 && conNum <= 3);
    assert(!isConUsed(conNum));
    cons[conNum].node = nodeIn;
    cons[conNum].flag = ConFlag::SET;
    cons[conNum].nodeNum = nodeIn->getNodeNum();
}

int BJT::genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    BJTModel *model = static_cast<BJTModel*>(getModel());
    Node *node = getNode(conNum);
    int B = getNodeNum(1) - 1, C = getNodeNum(0) - 1, E = getNodeNum(2) - 1, t = getNodeNum(conNum) - 1;
    double V_BE = 0, V_BC = 0;

    // 判断接地情况
    if (B != -1) {
        V_BE = (E != -1) ? (matrix->oldRHS[B] - matrix->oldRHS[E]) : matrix->oldRHS[B];
        V_BC = (C != -1) ? (matrix->oldRHS[B] - matrix->oldRHS[C]) : matrix->oldRHS[B];
    } else {
        V_BE = (E != -1) ? -matrix->oldRHS[E] : 0;
        V_BC = (C != -1) ? -matrix->oldRHS[C] : 0;
    }

    // 填充右端项
    switch (conNum) {
    case 0:
        matrix->RHS[t] += model->calculateIc(V_BE, V_BC);
        break;
    case 1:
        matrix->RHS[t] -= model->calculateIc(V_BE, V_BC) + model->calculateIe(V_BE, V_BC);
        break;
    case 2:
        matrix->RHS[t] += model->calculateIe(V_BE, V_BC);
        break;
    }
    return OK;
}

int BJT::genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    BJTModel *model = static_cast<BJTModel*>(getModel());
    // BJTType type = model->getType();
    int t = getNodeNum(conNum) - 1;
    int B = getNodeNum(1) - 1, C = getNodeNum(0) - 1, E = getNodeNum(2) - 1;
    double V_BE = 0, V_BC = 0;

    // 判断接地情况
    if (B != -1) {
        V_BE = (E != -1) ? (matrix->oldRHS[B] - matrix->oldRHS[E]) : matrix->oldRHS[B];
        V_BC = (C != -1) ? (matrix->oldRHS[B] - matrix->oldRHS[C]) : matrix->oldRHS[B];
    } else {
        V_BE = (E != -1) ? -matrix->oldRHS[E] : 0;
        V_BC = (C != -1) ? -matrix->oldRHS[C] : 0;
    }

    // 填充雅可比矩阵
    switch (conNum) {
    case 0: // C
        *matrix->elementMap[matrix->elementMapCountLocal++] += model->calculateCCJAC(V_BE, V_BC);
        // std::cout << model->calculateCCJAC(V_BE, V_BC) << std::endl;
        if (B != -1) *matrix->elementMap[matrix->elementMapCountLocal++] += model->calculateCBJAC(V_BE, V_BC);
        // std::cout << model->calculateCBJAC(V_BE, V_BC, B) << std::endl;
        if (E != -1) *matrix->elementMap[matrix->elementMapCountLocal++] += model->calculateCEJAC(V_BE, V_BC);
        // std::cout << model->calculateCEJAC(V_BE, V_BC, E) << std::endl;
        break;

    case 1: // B
        if (C != -1) *matrix->elementMap[matrix->elementMapCountLocal++] -= model->calculateBCJAC(V_BE, V_BC);
        // std::cout << model->calculateBCJAC(V_BE, V_BC, C) << std::endl;
        *matrix->elementMap[matrix->elementMapCountLocal++] -= model->calculateBBJAC(V_BE, V_BC);
        // std::cout << model->calculateBBJAC(V_BE, V_BC) << std::endl;
        if (E != -1) *matrix->elementMap[matrix->elementMapCountLocal++] -= model->calculateBEJAC(V_BE, V_BC);
        // std::cout << model->calculateBEJAC(V_BE, V_BC, E) << std::endl;
        break;
    
    case 2: // E
        if (C != -1) *matrix->elementMap[matrix->elementMapCountLocal++] += model->calculateECJAC(V_BE, V_BC);
        // std::cout << model->calculateECJAC(V_BE, V_BC, C) << std::endl;
        if (B != -1) *matrix->elementMap[matrix->elementMapCountLocal++] += model->calculateEBJAC(V_BE, V_BC);
        // std::cout << model->calculateEBJAC(V_BE, V_BC, B) << std::endl;
        *matrix->elementMap[matrix->elementMapCountLocal++] += model->calculateEEJAC(V_BE, V_BC);
        // std::cout << model->calculateEEJAC(V_BE, V_BC) << std::endl;
        break;
    }
    return OK;
}

int BJT::setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout,std::unique_ptr<Analysis>& analysis, const int conNum) {
    int t = getNodeNum(conNum) - 1, C = getNodeNum(0) - 1, B = getNodeNum(1) - 1, E = getNodeNum(2) - 1;

    // F(X)
    switch (conNum) {
    case 0:
        Fout[t] << " + " << getName() << "_Ic";
        break;
    case 1:
        Fout[t] << " - " << getName() << "_Ic"
                << " - " << getName() << "_Ie";
        break;
    case 2:
        Fout[t] << " + " << getName() << "_Ie";
        break;
    }

    // JAC
    BJTModel *model = static_cast<BJTModel*>(getModel());
    BJTType type = model->getType();
    switch (conNum) {
    case 0: // C
        if (type == BJTType::NPN) {
            // CC
            // - IS / alpha_r * v_T * (e ^ (V_BC * v_T))
            if (addCOOElement(coo, C, C, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - IS / alpha_r * v_T * exp((";
            if (B != -1) (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
            (*coo)->JACOutFile <<" - X_" << getNode(0)->getNodeName() << ") * v_T)";

            // CB
            // + IS / alpha_r * v_T * (e ^ (V_BC * v_T)) - IS * v_T * (e ^ (V_BE * v_T))
            if (B != -1) {
                if (addCOOElement(coo, C, B, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " + IS / alpha_r * v_T * exp(("
                                   << "X_" << getNode(1)->getNodeName()
                                   << " - X_" << getNode(0)->getNodeName() << ") * v_T)"
                                   << " - IS * v_T * exp((";
                (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
                if (E != -1) (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
            }
            
            // CE
            // + IS * v_T * (e ^ (V_BE * v_T)) 
            if (E != -1) {
                if (addCOOElement(coo, C, E, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " + IS * v_T * exp((";
                if (B != -1) (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
                (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName()
                                   << ") * v_T)";
            }

        } else {
            // CC
            if (addCOOElement(coo, C, C, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + IS / alpha_r * v_T * exp( - (";
            if (B != -1) (*coo)->JACOutFile<< "X_" << getNode(1)->getNodeName();
            (*coo)->JACOutFile <<" - X_" << getNode(0)->getNodeName() << ") * v_T)";

            // CB
            if (B != -1) {
                if (addCOOElement(coo, C, B, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - IS / alpha_r * v_T * exp( - ("
                                   << "X_" << getNode(1)->getNodeName()
                                   << " - X_" << getNode(0)->getNodeName() << ") * v_T)"
                                   << " + IS * v_T * exp( - (";
                (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
                if (E != -1) (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
            }

            // CE
            if (E != -1) {
                if (addCOOElement(coo, C, E, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - IS * v_T * exp( - (";
                if (B != -1) (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
                (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName()
                                   << ") * v_T)";
            }
        }
        break;
    
    case 1: // B
        if (type == BJTType::NPN) {
            // BC
            // + IS * v_T * (e ^ (V_BC * v_T)) - IS / alpha_r * v_T *(e ^ (V_BC * v_T))
            if (C != -1) {
                if (addCOOElement(coo, B, C, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - IS * v_T * exp((X_" << getNode(1)->getNodeName()
                                   << " - X_" << getNode(0)->getNodeName() << ") * v_T)"
                                   << " + IS / alpha_r * v_T * exp ((X_" << getNode(1)->getNodeName()
                                   << " - X_" << getNode(0)->getNodeName() << ") * v_T)";
            }

            // BB
            // + IS / alpha_f * v_T * (e ^ (V_BE * v_T)) - IS * v_T * (e ^ (V_BC * v_T)) +
            // IS / alpha_r * v_T * (e ^ (V_BC * v_T)) - IS * v_T * (e ^ (V_BE * v_T)) 
            if (addCOOElement(coo, B, B, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - IS / alpha_f * v_T * exp(("
                               << "X_" << getNode(1)->getNodeName();
            if (E != -1) (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";
            (*coo)->JACOutFile << " + IS * v_T * exp(("
                               << "X_" << getNode(1)->getNodeName();
            if (C != -1) (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";
            (*coo)->JACOutFile << " - IS / alpha_r * v_T * exp(("
                               << "X_" << getNode(1)->getNodeName();
            if (C != -1) (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";
            (*coo)->JACOutFile << " + IS * v_T * exp(("
                               << "X_" << getNode(1)->getNodeName();
            if (E != -1) (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";

            // BE
            // - IS / alpha_f * v_T * (e ^ (V_BE * v_T)) + IS * v_T * (e ^ (V_BE * v_T))
            if (E != -1) {
                if (addCOOElement(coo, B, E, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " + IS / alpha_f * v_T * exp(("
                        << "X_" << getNode(1)->getNodeName()
                        << " - X_" << getNode(2)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
                (*coo)->JACOutFile << " - IS * v_T * exp(("
                        << "X_" << getNode(1)->getNodeName()
                        << " - X_" << getNode(2)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
            }
            
        } else {
            // BC
            if (C != -1) {
                if (addCOOElement(coo, B, C, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " + IS * v_T * exp( - (X_" << getNode(1)->getNodeName();
                (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName() << ") * v_T)";
                (*coo)->JACOutFile << " - IS / alpha_r * v_T * exp ( - (X_" << getNode(1)->getNodeName();
                (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName() << ") * v_T)";
            }
            // BB
            if (addCOOElement(coo, B, B, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + IS / alpha_f * v_T * exp( - ("
                               << "X_" << getNode(1)->getNodeName();
            if (E != -1) (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";
            (*coo)->JACOutFile << " - IS * v_T * exp( - ("
                               << "X_" << getNode(1)->getNodeName();
            if (C != -1) (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";
            (*coo)->JACOutFile << " + IS / alpha_r * v_T * exp( - ("
                               << "X_" << getNode(1)->getNodeName();
            if (C != -1) (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";
            (*coo)->JACOutFile << " - IS * v_T * exp( - ("
                               << "X_" << getNode(1)->getNodeName();
            if (E != -1) (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";
            // BE
            if (E != -1) {
                if (addCOOElement(coo, B, E, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - IS / alpha_f * v_T * exp( - ("
                        << "X_" << getNode(1)->getNodeName()
                        << " - X_" << getNode(2)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
                (*coo)->JACOutFile << " + IS * v_T * exp( - ("
                        << "X_" << getNode(1)->getNodeName()
                        << " - X_" << getNode(2)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
            }
        }
        break;
    
    case 2: // E
        if (type == BJTType::NPN) {
            // EC
            // + IS * v_T * (e ^ (V_BC * v_T))
            if (C != -1) {
                if (addCOOElement(coo, E, C, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " + IS * v_T * exp((";
                if (B != -1) (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
                (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
            }

            // EB
            // + IS / alpha_f * v_T * (e ^ (V_BE * v_T)) - IS * v_T * (e ^ (V_BC * v_T))
            if (B != -1) {
                if (addCOOElement(coo, E, B, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " + IS / alpha_f * v_T * exp(("
                                   << "X_" << getNode(1)->getNodeName();
                if (E != -1) (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
                (*coo)->JACOutFile << " - IS * v_T * exp(("
                                   << "X_" << getNode(1)->getNodeName();
                if (C != -1) (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
            }

            // EE
            // - IS / alpha_f * v_T * (e ^ (V_BE * v_T))
            if (addCOOElement(coo, E, E, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - IS / alpha_f * v_T * exp((";
            if (B != -1) (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
            (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";

        } else {
            // EC
            if (C != -1) {
                if (addCOOElement(coo, E, C, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - IS * v_T * exp( - (";
                if (B != -1) (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
                (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
            }

            // EB
            if (B != -1) {
                if (addCOOElement(coo, E, B, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - IS / alpha_f * v_T * exp( - ("
                                   << "X_" << getNode(1)->getNodeName();
                if (E != -1) (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
                (*coo)->JACOutFile << " + IS * v_T * exp(-("
                                   << "X_" << getNode(1)->getNodeName();
                if (C != -1) (*coo)->JACOutFile << " - X_" << getNode(0)->getNodeName();
                (*coo)->JACOutFile << ") * v_T)";
            }

            // EE
            if (addCOOElement(coo, E, E, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + IS / alpha_f * v_T * exp( - (";
            if (B != -1) (*coo)->JACOutFile << "X_" << getNode(1)->getNodeName();
            (*coo)->JACOutFile << " - X_" << getNode(2)->getNodeName();
            (*coo)->JACOutFile << ") * v_T)";
        }
        break;
    }
    return OK;
}

void BJT::printMessage(std::ofstream &outFile, int conNum) {
    outFile << "        编号：" << getDeviceNum() << "    类型："
            << "BJT"
            << " 连接端口：" << conNum << "    名称：" << getName() << std::endl;
}

int Inductor::genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1;
    int id = analysis->getCirNodeTotal() + analysis->getCirBranchVsource() + getDeviceNum() - 1;

    switch (analysis->getCurMode()) {
    case AnalysisType::DC_Homotopy:
    case AnalysisType::DC_NR:
    case AnalysisType::DC_PTran:
        
        // // 判断接地情况
        // if (conNum == 0 && getNode(1)->getNodeNum() == 0 || conNum == 1 && getNode(0)->getNodeNum() == 0) {
        //     std::cerr << "Don't inductor link to ground!" << std::endl;
        //     return NA;
        // }

        if (conNum == 0) {
            // F_i = ... + X_id
            // F_id = X_i - X_j
            matrix->RHS[i] += matrix->oldRHS[id];
            if (j != -1)
                matrix->RHS[id] += matrix->oldRHS[i] - matrix->oldRHS[j];
            else
                matrix->RHS[id] += matrix->oldRHS[i];
        } else if (conNum == 1) {
            matrix->RHS[j] -= matrix->oldRHS[id];
        }
        break;
    case AnalysisType::TRAN:
        // F_i  = ... + x_id
        // F_id = X_id - X'_id - h / L * (V_i - V_j)
        if (conNum == 0) {
            matrix->RHS[i] += matrix->oldRHS[id];
            if (j != -1)
                matrix->RHS[id] += matrix->oldRHS[id] - matrix->oldoldRHS[id] - analysis->getTranStepSize() / getValue() * (matrix->oldRHS[i] - matrix->oldRHS[j]);
            else
                matrix->RHS[id] += matrix->oldRHS[id] - matrix->oldoldRHS[id] - analysis->getTranStepSize() / getValue() * matrix->oldRHS[i];
        } else if (conNum == 1) {
            matrix->RHS[j] -= matrix->oldRHS[id];
        }
        break;
    case AnalysisType::no_An:
        break;
    default:
        std::cerr << "Unknow Analysis Type!" << std::endl;
        return NA;
    }
    return OK;
}

int Inductor::genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int id = analysis->getCirNodeTotal() + analysis->getCirBranchVsource() + getDeviceNum() - 1;
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1;
    // 
    switch (analysis->getCurMode()) {
    case AnalysisType::DC_Homotopy:
    case AnalysisType::DC_NR:
    case AnalysisType::DC_PTran:
        if (conNum == 0) {
            // JAC(i, id) = 1
            *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            // JAC(id, i) = 1, JAC(id, j) = -1
            *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            if (j != -1)
                *matrix->elementMap[matrix->elementMapCountLocal++] -= 1;
        } else if (conNum == 1) {
            *matrix->elementMap[matrix->elementMapCountLocal++] -= 1;
        } else {
            std::cerr << "unkown conNum" << std::endl;
            return NA;
        }
        break;
    case AnalysisType::TRAN:
        // F_i  = ... + X_id
        // F_id = X_id - X'_id - h / L * (X_i - X_j)
        // JAC(i, id) += 1
        // JAC(id, id) += 1, JAC(id, i) -= h / L, JAC(id, j) += h / L;
        if (conNum == 0) {
            *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            *matrix->elementMap[matrix->elementMapCountLocal++] -= analysis->getTranStepSize() / getValue();
            if (j != -1)
                *matrix->elementMap[matrix->elementMapCountLocal++] += analysis->getTranStepSize() / getValue();
        } else if (conNum == 1) {
            // F_j = ... - X_id
            // JAC(j, id) -= 1;
            *matrix->elementMap[matrix->elementMapCountLocal++] -= 1;
        } else {
            std::cerr << "unkown conNum" << std::endl;
            return NA;
        }
        break;
    case AnalysisType::no_An:
        break;
    default:
        std::cerr << "Unknow Analysis Type!" << std::endl;
        return NA;
    }
    return OK;
}

int Inductor::setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    int id = analysis->getCirNodeTotal() + analysis->getCirBranchVsource() + getDeviceNum() - 1;

    // F(x)
    switch (analysis->getCurMode()) {
    case AnalysisType::DC_Homotopy:
    case AnalysisType::DC_NR:
    case AnalysisType::DC_PTran:
        if (conNum == 0) {
            Fout[t] << " + X_" << getName();
            Fout[id] << " + (X_" << getNode(0)->getNodeName();
            if (j != -1) Fout[id] << " - X_" << getNode(1)->getNodeName();
            Fout[id] << ")";
        } else if (conNum == 1) {
            Fout[t] << " - X_" << getName();
        } else {
            std::cerr << "unkown conNum" << std::endl;
            return NA;
        }
        break;
    case AnalysisType::TRAN:
        if (conNum == 0) {
            Fout[t] << " + X_" << getName();
            Fout[id] << " + X_" << getName() <<" - X'_" << getName() 
                     << " - h / L * (X_" << getNode(0)->getNodeName();
            if (j != -1) Fout[id] << " - X_" << getNode(1)->getNodeName();
            Fout[id] << ")";
        } else if (conNum == 1) {
            Fout[t] << " - X_" << getName();
        } else {
            std::cerr << "unkown conNum" << std::endl;
            return NA;
        }
        break;
    case AnalysisType::no_An:
        break;
    default:
        std::cerr << "unkown conNum" << std::endl;
        return NA;
    }

    // JAC
    switch (analysis->getCurMode()) {
    case AnalysisType::DC_Homotopy:
    case AnalysisType::DC_NR:
    case AnalysisType::DC_PTran:
        if (conNum == 0) {
            if (addCOOElement(coo, t, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + 1";
            if (addCOOElement(coo, id, i, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + 1";
            if (j != -1) {
                if (addCOOElement(coo, id, j, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - 1";
            }
        } else if (conNum == 1) {
            if (addCOOElement(coo, t, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - 1";
        } else {
            std::cerr << "unkown conNum" << std::endl;
            return NA;
        }
        break;
    case AnalysisType::TRAN:
        if (conNum == 0) {
            if (addCOOElement(coo, t, id, elementMapCount) == NA) return NA;;
            (*coo)->JACOutFile << " + 1";
            if (addCOOElement(coo, id, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + 1";
            if (addCOOElement(coo, id, i, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - h / L";
            if (j != -1) {
                if (addCOOElement(coo, id, j, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " + h / L";
            }
            
        } else if (conNum == 1) {
            if (addCOOElement(coo, t, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - 1";
        } else {
            std::cerr << "unkown conNum" << std::endl;
            return NA;
        }
        break;
    case AnalysisType::no_An:
        break;
    default:
        std::cerr << "Unknow Analysis Type!" << std::endl;
        return NA;
    }

    return OK;
}

void Inductor::printMessage(std::ofstream &outFile, int conNum) {
    outFile << "        编号：" << getDeviceNum() << "    类型："
            << "Inductor"
            << " 连接端口：" << conNum << "    名称：" << getName() << std::endl;
    outFile << "        value: " << getValue() << std::endl;
}

int Capacitor::genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    double V_ij_k_1 = 0, V_ij_k = 0;

    // 处理接地情况
    if (j == -1) {
        V_ij_k_1 = matrix->oldRHS[i];
        V_ij_k =  matrix->oldoldRHS[i];
    } else {
        V_ij_k_1 = matrix->oldRHS[i] - matrix->oldRHS[j];
        V_ij_k =  matrix->oldoldRHS[i] - matrix->oldoldRHS[j];
    }

    switch (analysis->getCurMode()) {
    case AnalysisType::no_An:
    case AnalysisType::DC_NR:
    case AnalysisType::DC_PTran:
    case AnalysisType::DC_Homotopy:
        break;
    case AnalysisType::TRAN:
        if (conNum == 0) {
            // F_i = ... + C / h * ((X_i - X_j) - (X'_i - X'_j))
            matrix->RHS[t] += getValue() / analysis->getTranStepSize() * (V_ij_k_1 - V_ij_k);
        } else {
            // F_j = ... - C / h * ((X_i - X_j) - (X'_i - X_j))
            matrix->RHS[t] -= getValue() / analysis->getTranStepSize() * (V_ij_k_1 - V_ij_k);
        }
        break;
    default:
        std::cerr << "Unknow Analysis Type!" << std::endl;
        return NA;
    }
    return OK;
}

int Capacitor::genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    double V_ij_k_1 = 0, V_ij_k = 0;

    // 处理接地情况
    if (j == -1) {
        V_ij_k_1 = matrix->oldRHS[i];
        V_ij_k =  matrix->oldoldRHS[i];
    } else {
        V_ij_k_1 = matrix->oldRHS[i] - matrix->oldRHS[j];
        V_ij_k =  matrix->oldoldRHS[i] - matrix->oldoldRHS[j];
    }

    switch (analysis->getCurMode()) {
    case AnalysisType::no_An:
    case AnalysisType::DC_NR:
    case AnalysisType::DC_PTran:
    case AnalysisType::DC_Homotopy:
        break;
    case AnalysisType::TRAN:
        if (conNum == 0) {
            // F_i = ... + C / h * ((X_i - X_j) - (X'_i - X'_j))
            *matrix->elementMap[matrix->elementMapCountLocal++] += getValue() / analysis->getTranStepSize();
            if (j != -1)
                *matrix->elementMap[matrix->elementMapCountLocal++] -= getValue() / analysis->getTranStepSize();
        } else {
            // F_j = ... - C / h * ((X_i - X_j) - (X'_i - X'_j))
            *matrix->elementMap[matrix->elementMapCountLocal++] -= getValue() / analysis->getTranStepSize();
            *matrix->elementMap[matrix->elementMapCountLocal++] += getValue() / analysis->getTranStepSize();
        }
        break;
    default:
        std::cerr << "Unknow Analysis Type!" << std::endl;
        return NA;
    }
    return OK;
}

int Capacitor::setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;

    // F(x)
    switch (analysis->getCurMode()) {
    case AnalysisType::no_An:
    case AnalysisType::DC_NR:
    case AnalysisType::DC_PTran:
    case AnalysisType::DC_Homotopy:
        break;
    case AnalysisType::TRAN:
        if (conNum == 0) {
            // F_i = ... + C / h * ((X_i - X_j) - (X'_i - X'_j))
            Fout[t] << " + " << getName() << " / h * ((X_" << getNode(0)->getNodeName();
            if (j != -1) Fout[t] << " - X_" << getNode(1)->getNodeName();
            Fout[t] << ") - (X'_" << getNode(0)->getNodeName();
            if (j != -1) Fout[t] << " - X'_" << getNode(1)->getNodeName();
            Fout[t] << "))";
        } else {
            // F_i = ... - C / h * ((X_i - X_j) - (X'_i - X'_j))
            Fout[t] << " - " << getName() << " / h * ((X_" << getNode(0)->getNodeName();
            Fout[t] << " - X_" << getNode(1)->getNodeName();
            Fout[t] << ") - (X'_" << getNode(0)->getNodeName();
            Fout[t] << " - X'_" << getNode(1)->getNodeName();
            Fout[t] << "))";
        }
        break;
    default:
        std::cerr << "Unknow Analysis Type!" << std::endl;
        return NA;
    }

    // JAC
    switch (analysis->getCurMode()) {
    case AnalysisType::no_An:
    case AnalysisType::DC_NR:
    case AnalysisType::DC_PTran:
    case AnalysisType::DC_Homotopy:
        break;
    case AnalysisType::TRAN:
        if (conNum == 0) {
            // F_i = ... + C / h * ((X_i - X_j) - (X'_i - X'_j))
            // JAC (t, i) += C / h, JAC(t, j) -= C / h
            if (addCOOElement(coo, t, i, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + C / h";
            if (j != -1) {
                if (addCOOElement(coo, t, j, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - C / h";
            }
        } else {
            // F_j = ... - C / h * ((X_i - X_j) - (X'_i - X'_j))
            if (addCOOElement(coo, t, i, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - C / h";
            if (addCOOElement(coo, t, j, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + C / h";
        }
        break;
    default:
        std::cerr << "Unknow Analysis Type!" << std::endl;
        return NA;
    }

    return OK;
}

void Capacitor::printMessage(std::ofstream &outFile, int conNum) {
    outFile << "        编号：" << getDeviceNum() << "    类型："
            << "Capacitor"
            << " 连接端口：" << conNum << "    名称：" << getName() << std::endl;
    outFile << "        value: " << getValue() << std::endl;
}

int Resistor::genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    switch (conNum) {
    case 0:
        matrix->RHS[t] += (j != -1) ? (matrix->oldRHS[i] - matrix->oldRHS[j]) / getValue() : matrix->oldRHS[i] / getValue();
        break;
    case 1:
        matrix->RHS[t] -= (matrix->oldRHS[i] - matrix->oldRHS[j]) / getValue();
        break;
    default:
        std::cerr << "Unknow ConNum!" << std::endl;
        return NA;
    }
    return OK;
}

int Resistor::genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    switch (conNum) {
    case 0:
        *matrix->elementMap[matrix->elementMapCountLocal++] += 1 / getValue();
        if (j != -1) *matrix->elementMap[matrix->elementMapCountLocal++] -= 1 / getValue();
        break;
    case 1:
        *matrix->elementMap[matrix->elementMapCountLocal++] -= 1 / getValue();
        *matrix->elementMap[matrix->elementMapCountLocal++] += 1 / getValue();
        break;
    default:
        std::cerr << "Unknow ConNum!" << std::endl;
        return NA;
    }
    return OK;
}

int Resistor::setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    
    // F(x)
    switch (conNum) {
    case 0:
        // F_i = ... + (X_i - X_j) / R
        Fout[t] << " + (X_" << getNode(0)->getNodeName();
        if (j != -1) Fout[t] << " - X_" << getNode(1)->getNodeName();
        Fout[t] << ") / " << getName();
        break;
    case 1:
        // F_j = ... - (X_i - X_j) / R
        Fout[t] << " - (X_" << getNode(0)->getNodeName();
        Fout[t] << " - X_" << getNode(1)->getNodeName();
        Fout[t] << ") / " << getName();
        break;
    default:
        std::cerr << "Unknow ConNum!" << std::endl;
        return NA;
    }

    // JAC
    switch (conNum) {
    case 0:
        // JAC(i, i) += 1 / R, JAC(i, j) -= 1 / R
        if (addCOOElement(coo, t, i, elementMapCount) == NA) return NA;
        (*coo)->JACOutFile << " + 1 / " << getName();
        if (j != -1) {
            if (addCOOElement(coo, t, j, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - 1 / " << getName();
        }
        break;
    case 1:
        // JAC(j, j) += 1 / R, JAC(j, i) -= 1 / R
        if (addCOOElement(coo, t, i, elementMapCount) == NA) return NA;
        (*coo)->JACOutFile << " - 1 / " << getName();
        if (addCOOElement(coo, t, j, elementMapCount) == NA) return NA;
        (*coo)->JACOutFile << " + 1 / " << getName();
        break;
    default:
        std::cerr << "Unknow ConNum!" << std::endl;
        return NA;
    }

    return OK;

}

void Resistor::printMessage(std::ofstream &outFile, int conNum) {
    outFile << "        编号：" << getDeviceNum() << "    类型："
            << "Resistor"
            << " 连接端口：" << conNum << "    名称：" << getName() << std::endl;
    outFile << "        value: " << getValue() << std::endl;
}

int Vsource::genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int id = analysis->getCirNodeTotal() + getDeviceNum() - 1;
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    double V = (j == -1) ? (matrix->oldRHS[i] - getValue()) : (matrix->oldRHS[i] - matrix->oldRHS[j] - getValue());

    if (analysis->getCurMode() == AnalysisType::DC_PTran) {
        switch (conNum) {
        case 0:
            // F_id = ... - X_id
            // DC: F_i = X_i - X_j - getValue()
            // X_i' = X_i - (X_id - X'_id) * L / h
            // DC_tran: F_i' = X_i - (X_id - X'_id) * L / h - X_j - getValue()
            // matrix->RHS[id] -= matrix->oldRHS[id];
            // matrix->RHS[t] += V - (matrix->oldRHS[id] - matrix->oldoldRHS[id]) * analysis->getPInductor() / analysis->getPTranStepSize();
            // F_i = ... + X_id
            // DC: F_id = X_i - X_j - getvalue()
            // DC_tran: F_id = X_i - (X_id - X'_id) * L / h - X_j - getValue()
            matrix->RHS[t] += matrix->oldRHS[id];
            matrix->RHS[id] += V - (matrix->oldRHS[id] - matrix->oldoldRHS[id]) * analysis->getPInductor() / analysis->getPTranStepSize();
            // matrix->RHS[id] += matrix->oldoldRHS[id] - matrix->oldRHS[id] + analysis->getPTranStepSize() / analysis->getPInductor() * V;
            break;
        case 1:
            // F_j = ... - X_id
            matrix->RHS[t] -= matrix->oldRHS[id];
            break;
        default:
            std::cerr << "Unknow ConNum!" << std::endl;
            return NA;
        }
    } else {
        switch (conNum) {
        case 0:
            // F_id = ... - X_id
            // F_i = X_i - X_j - getValue()
            // matrix->RHS[id] -= matrix->oldRHS[id];
            // matrix->RHS[t] += V;
            // F_i = ... + X_id
            // DC: F_id = X_i - X_j - getvalue()
            matrix->RHS[t] += matrix->oldRHS[id];
            matrix->RHS[id] += V;
            break;
        case 1:
            // F_j = ... - X_id
            matrix->RHS[t] -= matrix->oldRHS[id];
            break;
        default:
            std::cerr << "Unknow ConNum!" << std::endl;
            return NA;
        }
    }
    return OK;
}

int Vsource::genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    
    if (analysis->getCurMode() == AnalysisType::DC_PTran) {
        switch (conNum) {
        case 0:
            // JAC(i, id) += 1, JAC(id, i) += 1, JAC(id, j) -= 1, JAC(id, id) -= L / h
            *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            if (j != -1) *matrix->elementMap[matrix->elementMapCountLocal++] -= 1;
            *matrix->elementMap[matrix->elementMapCountLocal++] -= analysis->getPInductor() / analysis->getPTranStepSize();
            // JAC(i, id) += 1, JAC(id, i) += h / L, JAC(id, j) -= h / L, JAC(id, id) -= 1
            // *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            // *matrix->elementMap[matrix->elementMapCountLocal++] += analysis->getPTranStepSize() / analysis->getPInductor();
            // if (j != -1) *matrix->elementMap[matrix->elementMapCountLocal++] -= analysis->getPTranStepSize() / analysis->getPInductor();
            // *matrix->elementMap[matrix->elementMapCountLocal++] -= 1;
            break;
        case 1:
            // JAC(j, id) -= 1
            *matrix->elementMap[matrix->elementMapCountLocal++] -= 1;
            break;
        default:
            std::cerr << "Unknow ConNum!" << std::endl;
            return NA;
        }
    } else {
        switch (conNum) {
        case 0:
            // JAC(i, id) += 1, JAC(id, i) += 1, JAC(id, j) -= 1
            *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            *matrix->elementMap[matrix->elementMapCountLocal++] += 1;
            if (j != -1) *matrix->elementMap[matrix->elementMapCountLocal++] -= 1;
            break;
        case 1:
            // JAC(j, id) -= 1
            *matrix->elementMap[matrix->elementMapCountLocal++] -= 1;
            break;
        default:
            std::cerr << "Unknow ConNum!" << std::endl;
            return NA;
        }
    }
    return OK;
}

int Vsource::setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int id = analysis->getCirNodeTotal() + getDeviceNum() - 1;
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;

    // F(x)
    switch (conNum) {
    case 0:
        // F_id = ... - X_id
        // F_i = X_i - X_j - getValue()
        // F_i' = X_i - (X_id - X'_id) * L / h - X_j - getValue()
        // F_i = ... - X_id
        // DC: F_id = X_i - X_j - getvalue()
        // DC_tran: F_id = X_i - (X_id - X'_id) * L / h - X_j - getValue()
        // Fout[id] << " - X_v" << getName();
        // Fout[t] << " + X_" << getNode(0)->getNodeName();
        // if (j != -1) Fout[t] << " - X_" << getNode(1)->getNodeName();
        // Fout[t] << " - " << getName();
        // if (analysis->getCurMode() == AnalysisType::DC_PTran) {
        //     Fout[t] << " - (X_v" << getName() << " - X'_v" << getName() << ") * L / h";
        // }
        Fout[t] << " + X_" << getName();
        Fout[id] << " + X_" << getNode(0)->getNodeName();
        if (j != -1) Fout[id] << " - X_" << getNode(1)->getNodeName();
        Fout[id] << " - " << getName();
        if (analysis->getCurMode() == AnalysisType::DC_PTran) {
            Fout[id] << " - (X_" << getName() << " - X'_" << getName() << ") * L / h";
        }
        break;
    case 1:
        // F_j = ... - X_id
        Fout[t] << " - X_v" << getName();
        break;
    default:
        std::cerr << "Unknow ConNum!" << std::endl;
        return NA;
    }

    // JAC
    if (analysis->getCurMode() == AnalysisType::DC_PTran) {
        switch (conNum) {
        case 0:
            // JAC(id, id) -= 1
            // JAC(i, i) += 1, JAC(i, j) -= 1, JAC(i, id) -= L / h
            // JAC(i, id) += 1, JAC(id, i) += 1, JAC(id, j) -= 1, JAC(id, id) -= L / h
            if (addCOOElement(coo, i, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + 1";
            if (addCOOElement(coo, id, i, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + 1";
            if (j != -1) {
                if (addCOOElement(coo, id, j, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - 1";
            }
            if (addCOOElement(coo, id, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - L / h";

            // // JAC(i, id) += 1, JAC(id, i) += h / L, JAC(id, j) -= h / L, JAC(id, id) -= 1
            // if (addCOOElement(coo, i, id, elementMapCount) == NA) return NA;
            // (*coo)->JACOutFile << " + 1";
            // if (addCOOElement(coo, id, i, elementMapCount) == NA) return NA;
            // (*coo)->JACOutFile << " + h / L";
            // if (j != -1) {
            //     if (addCOOElement(coo, id, j, elementMapCount) == NA) return NA;
            //     (*coo)->JACOutFile << " - h / L";
            // }
            // if (addCOOElement(coo, id, id, elementMapCount) == NA) return NA;
            // (*coo)->JACOutFile << " - 1";
            break;
        case 1:
            // JAC(j, id) -= 1
            if (addCOOElement(coo, j, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - 1";
            break;
        default:
            std::cerr << "Unknow ConNum!" << std::endl;
            return NA;
        }
    } else {
        switch (conNum) {
        case 0:
            // JAC(id, id) -= 1
            // JAC(i, i) += 1, JAC(i, j) -= 1
            // JAC(i, id) += 1, JAC(id, i) += 1, JAC(id, j) -= 1
            if (addCOOElement(coo, i, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + 1";
            if (addCOOElement(coo, id, i, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " + 1";
            if (j != -1) {
                if (addCOOElement(coo, id, j, elementMapCount) == NA) return NA;
                (*coo)->JACOutFile << " - 1";
            }
            break;
        case 1:
            // JAC(j, id) -= 1
            if (addCOOElement(coo, j, id, elementMapCount) == NA) return NA;
            (*coo)->JACOutFile << " - 1";
            break;
        default:
            std::cerr << "Unknow ConNum!" << std::endl;
            return NA;
        }
    }
    return OK;

}

void Vsource::printMessage(std::ofstream &outFile, int conNum) {
    outFile << "        编号：" << getDeviceNum() << "    类型："
            << "VSource"
            << " 连接端口：" << conNum << "    名称：" << getName() << std::endl;
    outFile << "        value: " << getValue() << std::endl;
}

int Isource::genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;
    switch (conNum){
    case 0:
        // F_i = ... + I
        matrix->RHS[t] += getValue();
        break;
    case 1:
        // F_j = ... - I
        matrix->RHS[t] -= getValue();
        break;
    default:
        std::cerr << "Unknow ConNum!" << std::endl;
        return NA;
    }
    return OK;
}

int Isource::genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) {
    return OK;
}

int Isource::setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) {
    int i = getNodeNum(0) - 1, j = getNodeNum(1) - 1, t = getNodeNum(conNum) - 1;

    // F(X)
    switch (conNum){
    case 0:
        Fout[t] << " + I" << getName();
        break;
    case 1:
        // F_j = ... - I
        Fout[t] << " - I" << getName();
        break;
    default:
        std::cerr << "Unknow ConNum!" << std::endl;
        return NA;
    }
    return OK;
}

void Isource::printMessage(std::ofstream &outFile, int conNum) {
    outFile << "        编号：" << getDeviceNum() << "    类型："
            << "ISource"
            << " 连接端口：" << conNum << "    名称：" << getName() << std::endl;
    outFile << "        value: " << getValue() << std::endl;
}



/****************************************************************************************************/
/***********************************   DeviceHead 类定义       ***************************************/
/****************************************************************************************************/
/****************************************************************************************************/

DeviceHead::DeviceHead() {
    deviceList = nullptr;
    indCount = resCount = capCount = vsCount = isCount = qCount = 0;
}

DeviceHead::~DeviceHead() {
    Device *device = deviceList, *next;
    while (device != nullptr) {
        next = device->getNext();
        delete device;
        device = next;
    }
    deviceList = nullptr;
}

Device* DeviceHead::addDevice(DeviceType type, double valueIn, Model *modelIn, char *nameIn) {
    Device *device;
    switch (type) {
        case DeviceType::BJT:
            qCount++;
            device = new BJT(type, valueIn, qCount, modelIn, nameIn);
            deviceListByType[0].push_back(device);
            break;
        case DeviceType::VSource:
            vsCount++;
            device = new Vsource(type, valueIn, vsCount, modelIn, nameIn);
            deviceListByType[1].push_back(device);
            break;
        case DeviceType::Inductor:
            indCount++;
            device = new Inductor(type, valueIn, indCount, modelIn, nameIn);
            deviceListByType[2].push_back(device);
            break;
        case DeviceType::Resistor:
            resCount++;
            device = new Resistor(type, valueIn, resCount, modelIn, nameIn);
            deviceListByType[3].push_back(device);
            break;
        case DeviceType::Capacitor:
            capCount++;
            device = new Capacitor(type, valueIn, capCount, modelIn, nameIn);
            deviceListByType[4].push_back(device);
            break;
        case DeviceType::ISource:
            isCount++;
            device = new Isource(type, valueIn, isCount, modelIn, nameIn);
            deviceListByType[5].push_back(device);
            break;
        // default:
        //     perror("NA: unknow device type");
        //     exit(1);
        //     break;
    }  
    // 头插
    device->setNext(deviceList);
    deviceList = device;

    return deviceList;
}

int DeviceHead::getCount(DeviceType type) {
    switch (type) {
        case DeviceType::BJT:
            return qCount;
        case DeviceType::VSource:
            return vsCount;
        case DeviceType::Inductor:
            return indCount;
        case DeviceType::Resistor:
            return resCount;
        case DeviceType::Capacitor:
            return capCount;
        case DeviceType::ISource:
            return isCount;
        default:
            std::cerr << "NA: unknow device type!" << std::endl;
            return NA;
    }
    return OK;
}

std::vector<Device *>& DeviceHead::getDeviceListByType(DeviceType type) {
    switch (type) {
    case DeviceType::BJT:
        return deviceListByType[0];
    case DeviceType::VSource:
        return deviceListByType[1];
    case DeviceType::Inductor:
        return deviceListByType[2];
    case DeviceType::Resistor:
        return deviceListByType[3];
    case DeviceType::Capacitor:
        return deviceListByType[4];
    case DeviceType::ISource:
        return deviceListByType[5];
    }
}

// 获取对应编号的器件
Device* DeviceHead::getDevice(int deviceNum) {
    return deviceNum == 0? deviceList: [](Device *tmp0, int tmp1) -> Device* {
        for (Device *ptr = tmp0; ptr != nullptr; ptr = ptr->getNext()) {
            if (ptr->getDeviceNum() == tmp1) return ptr;
        }
        return nullptr;
    }(deviceList, deviceNum);
}


/****************************************************************************************************/
/*******************************       NodeHead 类定义        ***************************************/
/****************************************************************************************************/
/****************************************************************************************************/

void NodeHead::addGroundNode() {
    Node *tmp_node = new Node("0", 0, NodeType::general_Node);
    symtbl["0"] = 0;
    id2Name[0] = "0";
    tmp_node->setNext(nodeList);
    nodeList = tmp_node;
    return;
}

void NodeHead::addNode(char *nodeName, Device *device, int conNum, NodeType nodeType) {
    // if (symtbl[nodeName] != 0 || !strcmp(nodeName, "0")) return; // 节点已经存在
    if (nodeType == NodeType::general_Node) {
        Node *tmp_node;
        if (symtbl[nodeName] == 0 && strcmp(nodeName, "0")) {
            symtbl[nodeName] = ++nodeCount;
            id2Name[nodeCount] = nodeName;
            tmp_node = new Node(nodeName, nodeCount, nodeType);
            tmp_node->setNext(nodeList);
            nodeList = tmp_node;
        } else {
            tmp_node = getNode(nodeName);
        }

        //维护信息
        tmp_node->connect(conNum, device);
        device->connect(conNum, tmp_node);
        
    } else if (nodeType == NodeType::volatge_Node) {
        // 不需要新建节点，它是虚拟的 【这里conNum当作前置数量使用】
        int id = device->getDeviceNum() + conNum;
        setName2ID(device->getName(), id);
    } else if (nodeType == NodeType::current_Node) {
        // 不需要新建节点，它是虚拟的【这里conNum当作前置数量使用】
        int id = device->getDeviceNum() + conNum;
        setName2ID(device->getName(), id);
    }
}

int NodeHead::getCount() {
    return nodeCount;
}

// 获取指定节点编号的节点实例
Node* NodeHead::getNode(int nodeNum) {
    Node *ptr = nodeList;
    while (ptr != nullptr) {
        if (ptr->getNodeNum() == nodeNum) return ptr;
        ptr = ptr->getNext();
    }
    return nullptr;
}

// 获取指定节点名称的节点实例
Node* NodeHead::getNode(const char *nodeName) {
    return nodeName == nullptr? nodeList: [nodeName](Node *tmp_node) -> Node* {
        for (Node *ptr = tmp_node; ptr != nullptr; ptr = ptr->getNext()) {
            if (strcmp(ptr->getNodeName(), nodeName) == 0) return ptr;
        }
        return nullptr;
    }(nodeList);
}

// 获取头节点
Node* NodeHead::getNode() {
    return nodeList;
}

void NodeHead::printMessage(std::ofstream& outfile) {
    Node *ptr = nodeList;
    while (ptr != nullptr){
        outfile << "节点" << ptr->getNodeName() << ", 编号为：" << ptr->getNodeNum() 
                << ", 所连器件数为：" << ptr->getConsCount() << std::endl;
        Connections *cons_ptr = ptr->getConList();
        while (cons_ptr != nullptr) {
            cons_ptr->device->printMessage(outfile, cons_ptr->conNum);
            cons_ptr = cons_ptr->next;
        }
        ptr = ptr->getNext();
    }
}

