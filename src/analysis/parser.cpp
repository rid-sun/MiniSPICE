#include <iostream>
#include <fstream>
#include "minispice.hpp"

//前向声明
double stripString2Val(char *stringIn);
std::pair<std::string, double> stripString2IdVal(char *stringIn);


// 不支持windows行尾
int parser(std::unique_ptr<Netlist>& netlist, std::unique_ptr<Analysis>& analysis, std::string inFileName, std::string outFileName) {
    
    // 读写文件流
    std::ifstream inFile;
    std::ofstream outFile;

    // 分析过程中所需的变量
    char buf[BufLength], buf1[BufLength], buf2[BufLength], typebuf[BufLength], nameBuf[NameLength];
    char *bufPtr, *charPtr1, *charPtr2;
    char *conBuf1, *conBuf2, *conBuf3, *conBuf4;

    NodeHead& nodeHead = netlist->getNodeHead();
    ModelHead& modelHead = netlist->getModelHead();
    DeviceHead& deviceHead = netlist->getDeviceHead();


    // 0. 处理文件相关
    inFile.open(inFileName, std::ios::in);
    if (!inFile) {
        std::cerr << "Failed to open file for reading." << std::endl;
        return MINISPICEERROR;
    }

    // 1. 读取网表信息
    // 1.1 标题行
    inFile.getline(buf, BufLength);
    netlist->setTitle(buf);

    // 1.2 点语句识别
    while (inFile.getline(buf, BufLength)) {
        // std::cout << buf << std::endl;
        // if (!inFile.good()) {
        //     std::cerr << "An error occurred while reading the file." << std::endl;
        //     inFile.close();
        //     return MINISPICEERROR;
        // }
        if (inFile.gcount() >= 2 && buf[inFile.gcount() - 2] == '\r') {
            buf[inFile.gcount() - 2] = '\0';  // Remove the trailing '\r'
        }
        if (*buf == '\0') continue;
        strcpy(buf1, buf);
        strcpy(buf2, strtok(buf1, " "));
        
        // .model
        if (!strcmp(buf2, ".model")) {
            // 保存模型名
            if ((charPtr2 = strtok(nullptr, " ")) == nullptr) {
                std::cerr << "Unknown device name." << std::endl;
                inFile.close();
                return MINISPICEERROR;
            }
            strcpy(nameBuf, charPtr2);
            charPtr1 = strtok(nullptr, " ");
            if (charPtr1 == nullptr || strcmp(charPtr1, "PNP") && strcmp(charPtr1, "NPN")) {
                std::cerr << "Unknown Transistor Type." << std::endl;
                inFile.close();
                return MINISPICEERROR;
            }
            // 保存模型类型
            strcpy(typebuf, charPtr1);

            // 获取后续参数信息
            double alpha_f = NA, alpha_r = NA, v_T = NA, IS = NA;
            while ((charPtr1 = strtok(nullptr, " ")) != nullptr) {
                if ((charPtr1[0] == 'I') && (charPtr1[1] == 'S') && (charPtr1[2] == '=')) {
                    IS = stripString2Val(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'F') && (charPtr1[2] == '=')) {
                    alpha_f = stripString2Val(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'R') && (charPtr1[2] == '=')) {
                    alpha_r = stripString2Val(charPtr1);
                }
                if ((charPtr1[0] == 'T') && (charPtr1[1] == 'E') && (charPtr1[2] == '=')) {
                    v_T = stripString2Val(charPtr1);
                }
            }
            modelHead.getBJTModel(nameBuf, !strcmp(typebuf, "NPN") ? BJTType::NPN : BJTType::PNP, IS, alpha_f, alpha_r, v_T);
        
        } else if (!strcmp(buf2, ".op")) {
            if ((charPtr2 = strtok(nullptr, " ")) == nullptr) {
                std::cout << "---Default Op Type-NR." << std::endl;
                analysis->setAnalysisType(AnalysisType::DC_NR);
            } else {
                strcpy(typebuf, charPtr2);
                if (!strcmp(typebuf, "NR")) {
                    std::cout << "---Op Type-NR." << std::endl;
                    analysis->setAnalysisType(AnalysisType::DC_NR);
                } else if (!strcmp(typebuf, "ptran")) {
                    std::cout << "---Op Type-ptran." << std::endl;
                    analysis->setAnalysisType(AnalysisType::DC_PTran);
                    if ((charPtr2 = strtok(nullptr, " ")) != nullptr) {
                        analysis->setPtranCapactor(strtod(charPtr2, nullptr));
                        std::cout << "------pCapactor val " << analysis->getPCapactor() << std::endl;
                    }
                    if ((charPtr2 = strtok(nullptr, " ")) != nullptr) {
                        analysis->setPtranInductor(strtod(charPtr2, nullptr));
                        std::cout << "------pInductor val " << analysis->getPInductor() << std::endl;
                    }
                    if ((charPtr2 = strtok(nullptr, " ")) != nullptr) {
                        analysis->setPtranInitStepSize(strtod(charPtr2, nullptr));
                        std::cout << "------ptran_init_step val " << analysis->getInitPtranStepSize() << std::endl;
                    }
                } else if (!strcmp(typebuf, "homotopy")) {
                    analysis->setAnalysisType(AnalysisType::DC_Homotopy);
                } else {
                    std::cerr << "Unknown Op Type." << std::endl;
                    inFile.close();
                    return MINISPICEERROR;
                }
            }
        } else if (!strcmp(buf2, ".tran")) {
            analysis->setAnalysisType(AnalysisType::TRAN);
            charPtr1 = strtok(nullptr, " ");
            if (charPtr1 == nullptr) {
                std::cerr << "Unknown tran stop time." << std::endl;
                inFile.close();
                return MINISPICEERROR;
            }
            double tran_stop = strtod(charPtr1, nullptr);
            analysis->setTranStop(tran_stop);
        } else if (!strcmp(buf2, ".ic")) {
            while ((charPtr1 = strtok(nullptr, " ")) != nullptr) {
                std::pair<std::string, double> tmp = stripString2IdVal(charPtr1);
                netlist->setIC(tmp.first, tmp.second);
            }
        } else if (!strcmp(buf2, ".nodeset")) {
            while ((charPtr1 = strtok(nullptr, " ")) != nullptr) {
                std::pair<std::string, double> tmp = stripString2IdVal(charPtr1);
                netlist->setNodeSet(tmp.first, tmp.second);
            }
        } else if (*buf2 == '.') {
            std::cerr << "Unknown dot sentence: " << buf2 << std::endl;
            inFile.close();
            return MINISPICEERROR;
        }
    }
    inFile.close();
    inFile.clear();

    std::cout << "---Dot-sentence read over!" << std::endl;

    // 1.3 器件扫描
    inFile.open(inFileName, std::ios::in);
    if (!inFile) {
        std::cerr << "Failed to open file for reading." << std::endl;
        return MINISPICEERROR;
    }
    while (inFile.getline(buf, BufLength)) {
        // if (!inFile.good()) {
        //     std::cerr << "An error occurred while reading the file." << std::endl;
        //     inFile.close();
        //     return MINISPICEERROR;
        // }
        if (inFile.gcount() >= 2 && buf[inFile.gcount() - 2] == '\r') {
            buf[inFile.gcount() - 2] = '\0';  // Remove the trailing '\r'
        }
        if (*buf == '\0') continue;

        strcpy(buf1, buf);
        strcpy(buf2, strtok(buf1, " "));

        if (*buf2 == '*') continue;

        if (isalpha(*buf2)) {
            strcpy(nameBuf, buf2);
            if (strlen(nameBuf) > NameLength) {
                std::cerr << "Device name length greater than maximum length." << std::endl;
                inFile.close();
                return MINISPICEERROR;
            }

            Device *device_ptr;
            double val;
            Model *model_ptr;
            char *modelName;
            
            switch (*buf2) {
            case 'v':
            case 'V':
                conBuf1 = strtok(nullptr, " ");
                conBuf2 = strtok(nullptr, " ");
                val = strtod(strtok(nullptr, " "), nullptr);
                device_ptr = deviceHead.addDevice(DeviceType::VSource, val, nullptr, nameBuf);
                if (!strcmp(conBuf1, "0")) {
                    std::cerr << "Positive can't connect to ground." << std::endl;
                    inFile.close();
                    return MINISPICEERROR;
                }
                nodeHead.addNode(conBuf1, device_ptr, 0, NodeType::general_Node);
                nodeHead.addNode(conBuf2, device_ptr, 1, NodeType::general_Node);
                break;

            case 'i':
            case 'I':
                conBuf1 = strtok(nullptr, " ");
                conBuf2 = strtok(nullptr, " ");
                val = strtod(strtok(nullptr, " "), nullptr);
                device_ptr = deviceHead.addDevice(DeviceType::ISource, val, nullptr, nameBuf);
                if (!strcmp(conBuf1, "0")) {
                    std::cerr << "Positive can't connect to ground." << std::endl;
                    inFile.close();
                    return MINISPICEERROR;
                }
                nodeHead.addNode(conBuf1, device_ptr, 0, NodeType::general_Node);
                nodeHead.addNode(conBuf2, device_ptr, 1, NodeType::general_Node);
                break;

            case 'q':
            case 'Q':
                conBuf1 = strtok(nullptr, " ");
                conBuf2 = strtok(nullptr, " ");
                conBuf3 = strtok(nullptr, " ");
                modelName = strtok(nullptr, " ");
                model_ptr = modelHead.getModel(modelName);
                device_ptr = deviceHead.addDevice(DeviceType::BJT, NA, model_ptr, nameBuf);
                if (!strcmp(conBuf2, "0")) {
                    std::cerr << "Base can't connect to ground." << std::endl;
                    inFile.close();
                    return MINISPICEERROR;
                }
                nodeHead.addNode(conBuf1, device_ptr, 0, NodeType::general_Node);
                nodeHead.addNode(conBuf2, device_ptr, 1, NodeType::general_Node);
                nodeHead.addNode(conBuf3, device_ptr, 2, NodeType::general_Node);
                break;

            case 'r':
            case 'R':
                conBuf1 = strtok(nullptr, " ");
                conBuf2 = strtok(nullptr, " ");
                val = strtod(strtok(nullptr, " "), nullptr);
                device_ptr = deviceHead.addDevice(DeviceType::Resistor, val, nullptr, nameBuf);
                if (!strcmp(conBuf1, "0")) {
                    std::cerr << "Positive can't connect to ground." << std::endl;
                    inFile.close();
                    return MINISPICEERROR;
                }
                nodeHead.addNode(conBuf1, device_ptr, 0, NodeType::general_Node);
                nodeHead.addNode(conBuf2, device_ptr, 1, NodeType::general_Node);
                break;

            case 'c':
            case 'C':
                conBuf1 = strtok(nullptr, " ");
                conBuf2 = strtok(nullptr, " ");
                val = strtod(strtok(nullptr, " "), nullptr);
                device_ptr = deviceHead.addDevice(DeviceType::Capacitor, val, nullptr, nameBuf);
                if (!strcmp(conBuf1, "0")) {
                    std::cerr << "Positive can't connect to ground." << std::endl;
                    inFile.close();
                    return MINISPICEERROR;
                }
                nodeHead.addNode(conBuf1, device_ptr, 0, NodeType::general_Node);
                nodeHead.addNode(conBuf2, device_ptr, 1, NodeType::general_Node);
                break;

            case 'l':
            case 'L':
                conBuf1 = strtok(nullptr, " ");
                conBuf2 = strtok(nullptr, " ");
                val = strtod(strtok(nullptr, " "), nullptr);
                device_ptr = deviceHead.addDevice(DeviceType::Inductor, val, nullptr, nameBuf);
                if (!strcmp(conBuf1, "0")) {
                    std::cerr << "Positive can't connect to ground." << std::endl;
                    inFile.close();
                    return MINISPICEERROR;
                }
                nodeHead.addNode(conBuf1, device_ptr, 0, NodeType::general_Node);
                nodeHead.addNode(conBuf2, device_ptr, 1, NodeType::general_Node);
                break;
            
            default:
                std::cerr << "Unknown device type." << std::endl;
                inFile.close();
                return MINISPICEERROR;
            }
        }
    }
    inFile.close();
    inFile.clear();

    std::cout << "---Device scanning over!" << std::endl;

    // 2. 处理信息
    analysis->setCirNodeTotal(nodeHead.getCount());
    analysis->setCirBranchVsource(deviceHead.getCount(DeviceType::VSource));
    analysis->setCirBranchInductor(deviceHead.getCount(DeviceType::Inductor));

    std::cout << "---Device connecting over!" << std::endl;

    // 3. 处理虚拟节点
    std::vector<Device *>& vscVector = deviceHead.getDeviceListByType(DeviceType::VSource);
    for (int i = 0; i < vscVector.size(); i++) {
        nodeHead.addNode(nullptr, vscVector[i], 0, NodeType::volatge_Node);
    }
    std::vector<Device *>& indVector = deviceHead.getDeviceListByType(DeviceType::Inductor);
    for (int i = 0; i < indVector.size(); i++) {
        nodeHead.addNode(nullptr, indVector[i], vscVector.size(), NodeType::current_Node);
    }
    if (vscVector.size() != deviceHead.getCount(DeviceType::VSource) || indVector.size() != deviceHead.getCount(DeviceType::Inductor)) {
        std::cerr << "Num of vsource or inductor don't  be equal to deviceHead." << std::endl;
        return MINISPICEERROR;
    }

    std::cout << "---Device virtual node setup over!" << std::endl;

    // 4. 打印输出信息
    outFile.open(outFileName + "_parser.txt", std::ios::out);
    if (!outFile) {
        std::cerr << "Failed to open file for writing." << std::endl;
        return MINISPICEERROR;
    }
    outFile << "Title: " << netlist->getTitle() << std::endl;
    outFile << "datum = " << netlist->getDatum() << "     cktNodeTotal = " << nodeHead.getCount() << std::endl;
    nodeHead.printMessage(outFile);

    std::cout << "---Message outputing over!" << std::endl;

    outFile.close();
    outFile.clear();

    return MINISPICEOK;
}


// 将传入的字符串 params=val 解析为 val
double stripString2Val(char *stringIn) {
    char buf[BufLength], buf2[BufLength];
    int a, b;
    strcpy(buf, stringIn);
    for (a = 0; buf[a] != '='; a++);
    a++;
    for (b = 0; buf[a] != '\0'; b++, a++)
        buf2[b] = buf[a];
    buf2[b] = '\0';
    return atof(buf2); //转成浮点数
}

// 将传入的字符串 v(id)=val 解析为 <id, val>
std::pair<std::string, double> stripString2IdVal(char *stringIn) {
    char buf[BufLength], buf2[BufLength], buf3[BufLength];
    int a, b;
    strcpy(buf, stringIn);
    for (a = 0; buf[a] != '('; a++);
    a++;
    for (b = 0; buf[a] != ')'; b++, a++)
        buf2[b] = buf[a];
    buf2[b] = '\0';
    for (; buf[a] != '='; a++);
    a++;
    for (b = 0; buf[a] != '\0'; b++, a++)
        buf3[b] = buf[a];
    buf3[b] = '\0';
    return std::pair<std::string, double>(std::string(buf2), strtod(buf3, nullptr));
}
