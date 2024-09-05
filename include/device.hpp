
#ifndef MINISPICEDEVICE_H
#define MINISPICEDEVICE_H

#include <memory>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <vector>
#include <unordered_map>
#include "MiniSPICEMatrix.hpp"

/* 这里主要定义了：五种枚举类型{器件类型、BJT类型、端口状态标志、逻辑类型、分析类型}
                  两种信息类结构体{器件端口信息、连接关系信息}
                  六种封装类{器件类、器件链表类、节点类、节点链表类、模型类、网表类}
*/

// 常量声明
constexpr int NameLength = 80;
extern const int BufLength;
extern const int NA;
extern const int OK;
// const int NameLength;
// const int BufLength;
// const int NA;
// const int OK;

// 器件类型
enum class DeviceType {
    BJT,
    VSource,
    ISource,
    Inductor,
    Resistor,
    Capacitor
}; 

// 晶体管类型
enum class BJTType {
    NPN,
    PNP
};

// 器件端口连接状态标志
enum class ConFlag {
    UNSET,
    SET
};

// 非线性电压/电流分支状态标志
// 避免重复添加伪元素器件
enum class PtranFlag {
    UNSET,
    SET
};

// 分析类型
enum class AnalysisType {
    no_An,
    DC_NR,
    DC_Homotopy,
    DC_PTran,
    TRAN
};

// 节点类型
enum class NodeType {
    general_Node,
    volatge_Node,
    current_Node
};

// 七种封装类
class Device;
class DeviceHead;
class Node;
class NodeHead;
class Model;
class BJTModel;
class ModelHead;
class Netlist;
class Analysis;

// 器件端口信息【器件->节点关系】
struct Connectors {
    ConFlag flag;             // 端口状态
    Node *node;               // 端口所连节点
    int nodeNum;              // 端口所连节点编号
    Connectors():flag(ConFlag::UNSET), node(nullptr), nodeNum(NA) {}
};

// 连接关系信息【节点->器件关系】
struct Connections {
    Connections *next;        // 自身预留接口，以成链
    Device *device;           // 指向所连器件
    int conNum;               // 该关系所连器件端口号
    Connections():next(nullptr), device(nullptr), conNum(NA) {}
};


/********************************************************************************************************/
/*******************************************DEVICE BEGIN*************************************************/

// 器件类抽象基类
class Device {
private:
    Connectors cons[4];       // 器件的4个端口
    Model *model;             // 器件对应模型
    DeviceType type;          // 器件类型
    Device *next;             // 自身预留接口，以成链
    int deviceNum;            // 器件实例在当前类型中的编号
    double value;             // 器件送值
    char name[NameLength];    // 器件名称

public:
    Device(DeviceType typeIn, double valueIn, int deviceNumIn, Model *modelIn, char *nameIn) {
        type = typeIn;
        value = valueIn;
        deviceNum = deviceNumIn;
        model = modelIn;
        strcpy(name, nameIn);
    }

    virtual ~Device() = default;

    /* 用来构建相应的方程 */
    virtual int genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) = 0;
    virtual int genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) = 0;
    virtual int setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) = 0;


    /* 打印器件信息 */
    virtual void printMessage(std::ofstream &outFile, int conNum) = 0;

    /* getter */
    DeviceType getType() {
        return type;
    }

    Device *getNext() {
        return next;
    }

    Node *getNode(int conNum) {
        return cons[conNum].node;
    }

    int getDeviceNum() {
        return deviceNum;
    }

    char *getName() {
        return name;
    }

    int getNodeNum(int conNum) {
        return cons[conNum].nodeNum;
    }

    double getValue() {
        return value;
    }

    Model* getModel() {
        return model;
    }

    /* setter */
    void setNext(Device *nextIn) {
        next = nextIn;
    }

    // 将器件端口连接到对应节点
    void connect(int conNum, Node *nodeIn);

    // 判断当前器件端口是否可用
    bool isConUsed(int conNum) {
        return cons[conNum].flag == ConFlag::SET ? true : false;
    }

    // // 在处理雅可比矩阵时添加临时元素
    // int addCOOElement(COOElement **head, int row, int col, int &elementMapCount) {
    //     int id = elementMapCount++;
    //     COOElement *ptr = new COOElement(id, row, col);
    //     if (ptr == nullptr) return NA;
    //     ptr->next = *head;
    //     *head = ptr;
    //     return 0;
    // }

};

// BJT 【con0：集电极C，con1：基极B，con2：发射级E】
class BJT: public Device {
public:
    BJT(DeviceType typeIn, double valueIn, int deviceNum, Model *modelIn, char *nameIn): 
        Device(typeIn, valueIn, deviceNum, modelIn, nameIn) {}

    ~BJT() {};

    /* 生成雅可比矩阵即KCL方程 */
    int genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout,std::unique_ptr<Analysis>& analysis, const int conNum) override;

    void printMessage(std::ofstream &outFile, int conNum) override;
    
};

// Inductor
class Inductor: public Device {
public:
    Inductor(DeviceType typeIn, double valueIn, int deviceNum, Model *modelIn, char *nameIn): 
        Device(typeIn, valueIn, deviceNum, modelIn, nameIn){}

    ~Inductor() {};

    /* 生成雅可比矩阵即KCL方程 */
    int genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    void printMessage(std::ofstream &outFile, int conNum) override;

};

// Capacitor
class Capacitor: public Device {
public:
    Capacitor(DeviceType typeIn, double valueIn, int deviceNum, Model *modelIn, char *nameIn): 
        Device(typeIn, valueIn, deviceNum, modelIn, nameIn){}

    ~Capacitor() {};

    /* 生成雅可比矩阵即KCL方程 */
    int genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    void printMessage(std::ofstream &outFile, int conNum) override;

};

// Resistor
class Resistor: public Device {
public:
    Resistor(DeviceType typeIn, double valueIn, int deviceNum, Model *modelIn, char *nameIn): 
        Device(typeIn, valueIn, deviceNum, modelIn, nameIn){}
    ~Resistor() {};

    /* 生成雅可比矩阵即KCL方程 */
    int genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    void printMessage(std::ofstream &outFile, int conNum) override;

};

// Vsource
class Vsource: public Device {
public:
    Vsource(DeviceType typeIn, double valueIn, int deviceNum, Model *modelIn, char *nameIn): 
        Device(typeIn, valueIn, deviceNum, modelIn, nameIn){}
    ~Vsource() {};

    /* 生成雅可比矩阵即KCL方程 */
    int genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    void printMessage(std::ofstream &outFile, int conNum) override;

};

// Isource
class Isource: public Device {
public:
    Isource(DeviceType typeIn, double valueIn, int deviceNum, Model *modelIn, char *nameIn): 
        Device(typeIn, valueIn, deviceNum, modelIn, nameIn){}

    ~Isource() {};

    /* 生成雅可比矩阵即KCL方程 */
    int genKCLEquation(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    int genKCLJAC(std::unique_ptr<SMPMatrix>& matrix, std::unique_ptr<Analysis>& analysis, const int conNum) override;
    
    int setup(COOElement **coo, int &elementMapCount, std::vector<std::ostringstream>& Fout, std::unique_ptr<Analysis>& analysis, const int conNum) override;

    void printMessage(std::ofstream &outFile, int conNum) override;

};


// 器件工厂
class DeviceHead {
private:
    DeviceHead();

    ~DeviceHead();

    // static DeviceHead deviceHead;
    Device *deviceList; // 器件链表
    std::vector<Device *> deviceListByType[6];
    int indCount, resCount, capCount, vsCount, isCount, qCount; // 每种器件的总数目

public:
    Device* addDevice(DeviceType type, double valueIn, Model *modelIn, char *nameIn);


    int getCount(DeviceType type);

    std::vector<Device *>& getDeviceListByType(DeviceType type);

    // 获取对应编号的器件
    Device *getDevice(int deviceNum=0);

    // 删除拷贝构造函数和赋值运算符确保单例实例的唯一性
    DeviceHead(const DeviceHead&) = delete;
    DeviceHead& operator=(const DeviceHead&) = delete;
    
    static DeviceHead& getInstance() {
        static DeviceHead deviceHead;
        return deviceHead;
    }

};

/*******************************************DEVICE  END**************************************************/
/********************************************************************************************************/



/********************************************************************************************************/
/*******************************************MODEL BEGIN**************************************************/

// 模型类抽象基类
class Model {
private:
    Model *next;
    char name[NameLength];             // 模型名称

public:
    Model(char *nameIn) {
        if (strlen(nameIn) > NameLength) {
            perror("Error, BJTModel initialize, exceed maxlength");
            exit(1);
        }
        strcpy(name, nameIn);
        next = nullptr;
    }

    virtual ~Model() = default;

    Model *getNext() {
        return next;
    }

    void setNext(Model *nextIn) {
        next = nextIn;
    }

    char *getName() {
        return name;
    }

};

// Ebers-Moll model
class BJTModel: public Model {
private:
    double IS, alpha_f, alpha_r, v_T;  // 模型参数
    BJTType type;                      // 模型类型

public:
    BJTModel(char *nameIn, BJTType typeIn, double ISin, double alpha_f_in, double alpha_r_in, double v_T_in): Model(nameIn) {
        type = typeIn;
        IS = ISin == NA ? -1e-16 : ISin;
        alpha_f = alpha_f_in == NA ? 0.99 : alpha_f_in;
        alpha_r = alpha_r_in == NA ? 0.5 : alpha_r_in;
        v_T = v_T_in == NA ? 38.78 : v_T_in;
    }

    ~BJTModel() override {}

    // 特性方程相关
    double calculateIe(double V_BE, double V_BC) {
        /* 这里的IE方向流向基区 */
        // NPN：IE = IS / alpha_f * ( e ^ ( V_BE * v_T) - 1 ) - IS * ( e ^ ( V_BC * v_T) - 1 )
        // PNP: IE = IS / alpha_f * ( e ^ ( V_EB * v_T) - 1 ) - IS * ( e ^ ( V_CB * v_T) - 1 )
        return type == BJTType::NPN ? (IS / alpha_f * (exp(V_BE * v_T) - 1) - IS * (exp(V_BC * v_T) - 1)) : (IS / alpha_f * (exp(-V_BE * v_T) - 1) - IS * (exp(-V_BC * v_T) - 1));
    }

    double calculateIc(double V_BE, double V_BC) {
        /* 这里的IC方向流向基区 */
        // NPN：IC = IS / alpha_r * ( e ^ ( V_BC * v_T) - 1 ) - IS * ( e ^ ( V_BE * v_T) - 1 )
        // PNP: IC = IS / alpha_r * ( e ^ ( V_CB * v_T) - 1 ) - IS * ( e ^ ( V_EB * v_T) - 1 ) 
        return type == BJTType::NPN ? (IS / alpha_r * (exp(V_BC * v_T) - 1) - IS * (exp(V_BE * v_T) - 1)) : (IS / alpha_r * (exp(-V_BC * v_T) - 1) - IS * (exp(-V_BE * v_T) - 1));
    }

    double calculateCEJAC(double V_BE, double V_BC) {
        // IS / alpha_r * ( e ^ ( V_BC * v_T) - 1 ) - IS * ( e ^ ( V_BE * v_T) - 1 )
        // + IS * v_T * (e ^ (V_BE * v_T)) 
        if (type == BJTType::NPN) {
            return IS * v_T * exp(V_BE * v_T);
        } else {
            return -IS * v_T * exp(-V_BE * v_T);
        }
    }

    double calculateCBJAC(double V_BE, double V_BC) {
        // IS / alpha_r * ( e ^ ( V_BC * v_T) - 1 ) - IS * ( e ^ ( V_BE * v_T) - 1 )
        // + IS / alpha_r * v_T * (e ^ (V_BC * v_T)) - IS * v_T * (e ^ (V_BE * v_T))
        if (type == BJTType::NPN) {
            return IS / alpha_r * v_T * exp(V_BC * v_T) - IS * v_T * exp(V_BE * v_T);
        } else {
            return -IS / alpha_r * v_T * exp(-V_BC * v_T) + IS * v_T * exp(-V_BE * v_T);
        }
    }

    double calculateCCJAC(double V_BE, double V_BC) {
        // IS / alpha_r * ( e ^ ( V_BC * v_T) - 1 ) - IS * ( e ^ ( V_BE * v_T) - 1 )
        // - IS / alpha_r * v_T * (e ^ (V_BC * v_T))
        if (type == BJTType::NPN) {
            return -IS / alpha_r * v_T * exp(V_BC * v_T);
        } else {
            return IS / alpha_r * v_T * exp(-V_BC * v_T);
        }
    }

    double calculateBEJAC(double V_BE, double V_BC) {
        // IS / alpha_f * ( e ^ ( V_BE * v_T) - 1 ) - IS * ( e ^ ( V_BC * v_T) - 1 ) +
        // IS / alpha_r * ( e ^ ( V_BC * v_T) - 1 ) - IS * ( e ^ ( V_BE * v_T) - 1 )
        // - IS / alpha_f * v_T * (e ^ (V_BE * v_T)) + IS * v_T * (e ^ (V_BE * v_T))
        if (type == BJTType::NPN) {
            return -IS / alpha_f * v_T * exp(V_BE * v_T) + IS * v_T * exp(V_BE * v_T);
        } else {
            return IS / alpha_f * v_T * exp(-V_BE * v_T) - IS * v_T * exp(-V_BE * v_T);
        }
    }

    double calculateBBJAC(double V_BE, double V_BC) {
        // IS / alpha_f * ( e ^ ( V_BE * v_T) - 1 ) - IS * ( e ^ ( V_BC * v_T) - 1 ) +
        // IS / alpha_r * ( e ^ ( V_BC * v_T) - 1 ) - IS * ( e ^ ( V_BE * v_T) - 1 )
        // + IS / alpha_f * v_T * (e ^ (V_BE * v_T)) - IS * v_T * (e ^ (V_BC * v_T)) +
        // IS / alpha_r * v_T * (e ^ (V_BC * v_T)) - IS * v_T * (e ^ (V_BE * v_T)) 
        if (type == BJTType::NPN) {
            return IS / alpha_f * v_T * exp(V_BE * v_T) - IS * v_T * exp(V_BC * v_T) + 
                   IS / alpha_r * v_T * exp(V_BC * v_T) - IS * v_T * exp(V_BE * v_T);
        } else {
            return -IS / alpha_f * v_T * exp(-V_BE * v_T) + IS * v_T * exp(-V_BC * v_T) - 
                   IS / alpha_r * v_T * exp(-V_BC * v_T) + IS * v_T * exp(-V_BE * v_T);
        }
    }

    double calculateBCJAC(double V_BE, double V_BC) {
        // IS / alpha_f * ( e ^ ( V_BE * v_T) - 1 ) - IS * ( e ^ ( V_BC * v_T) - 1 ) +
        // IS / alpha_r * ( e ^ ( V_BC * v_T) - 1 ) - IS * ( e ^ ( V_BE * v_T) - 1 )
        // + IS * v_T * (e ^ (V_BC * v_T)) - IS / alpha_r * v_T *(e ^ (V_BC * v_T))
        if (type == BJTType::NPN) {
            return IS * v_T * exp(V_BC * v_T) - IS / alpha_r * v_T * exp(V_BC * v_T);
        } else {
            return -IS * v_T * exp(-V_BC * v_T) + IS / alpha_r * v_T * exp(-V_BC * v_T);
        }
    }

    double calculateECJAC(double V_BE, double V_BC) {
        // IS / alpha_f * ( e ^ ( V_BE * v_T) - 1 ) - IS * ( e ^ ( V_BC * v_T) - 1 )
        // + IS * v_T * (e ^ (V_BC * v_T))
        if (type == BJTType::NPN) {
            return IS * v_T * exp(V_BC * v_T);
        } else {
            return -IS * v_T * exp(-V_BC * v_T);
        }
    }

    double calculateEBJAC(double V_BE, double V_BC) {
        // IS / alpha_f * ( e ^ ( V_BE * v_T) - 1 ) - IS * ( e ^ ( V_BC * v_T) - 1 )
        // + IS / alpha_f * v_T * (e ^ (V_BE * v_T)) - IS * v_T * (e ^ (V_BC * v_T))
        if (type == BJTType::NPN) {
            return IS / alpha_f * v_T * exp(V_BE * v_T) - IS * v_T * exp(V_BC * v_T);
        } else {
            return -IS / alpha_f * v_T * exp(-V_BE * v_T) + IS * v_T * exp(-V_BC * v_T);
        }
    }

    double calculateEEJAC(double V_BE, double V_BC) {
        // IS / alpha_f * ( e ^ ( V_BE * v_T) - 1 ) - IS * ( e ^ ( V_BC * v_T) - 1 )
        // - IS / alpha_f * v_T * (e ^ (V_BE * v_T))
        if (type == BJTType::NPN) {
            return -IS / alpha_f * v_T * exp(V_BE * v_T);
        } else {
            return IS / alpha_f * v_T * exp(-V_BE * v_T);
        }
    }


    BJTType getType() {
        return type;
    }

};


// 模型工厂
class ModelHead {
private:
    ModelHead() {
        modelList = nullptr;
        modelCount = 0;
    }

    ~ModelHead() {
        Model *model = modelList, *next;
        while (model != nullptr) {
            next = model->getNext();
            delete model;
            model = next;
        }
        modelList = nullptr;
    }

    // static ModelHead modelHead;
    Model *modelList;            // 器件链表
    int modelCount;              // 每种器件的总数目

public:
    Model *getBJTModel(char *nameIn, BJTType typeIn, double ISin, double alpha_f_in, double alpha_r_in, double v_T_in) {
        Model *ptr = new BJTModel(nameIn, typeIn, ISin, alpha_f_in, alpha_r_in, v_T_in);
        ptr->setNext(modelList);
        modelList = ptr;
        modelCount++;
        return modelList;
    }

    int getCount() {
       return modelCount;
    }

    /* 匿名表达式 */
    Model *getModel(char *modelName = nullptr) {
        return modelName == nullptr? modelList: [modelName](Model *tmp_model) -> Model* {
            for (Model *ptr = tmp_model; ptr != nullptr; ptr = ptr->getNext()) {
                if (strcmp(ptr->getName(), modelName) == 0) return ptr;
            }
            return nullptr;
        }(modelList);
    }

    // 删除拷贝构造函数和赋值运算符确保单例实例的唯一性
    ModelHead(const ModelHead&) = delete;
    ModelHead& operator=(const ModelHead&) = delete;
    
    static ModelHead& getInstance() {
        static ModelHead modelHead;
        return modelHead;
    }

};


/*******************************************MODELL  END**************************************************/
/********************************************************************************************************/


/********************************************************************************************************/
/*******************************************NODE  BEGIN**************************************************/

// 电路节点类
class Node {
private:
    Node *next;                   // 自身预留接口，以成链
    Connections *conList;         // 节点所相关的连接关系
    char nodeName[NameLength];    // 节点名称
    int nodeNum;                  // 节点编号
    // PtranFlag ptranFlag;          // Ptran Flag
    int consCount;                // 连接关系数目
    NodeType nodeType;            // 节点类型

public:
    Node(const char *name, int num, NodeType nodeTypeIn) {
        if (strlen(name) > NameLength) {
            perror("Error, BJTModel initialize, exceed maxlength");
            exit(1);
        }
        strcpy(nodeName, name);
        nodeNum = num;
        next = nullptr;
        conList = nullptr;
        // ptranFlag = PtranFlag::UNSET;
        consCount = 0;
        nodeType = nodeTypeIn;
    }

    ~Node() {
        Connections *ptr = conList, *next;
        while (ptr != nullptr) {
            next = ptr->next;
            delete ptr;
            ptr = next;
        }
        conList = nullptr;
    };

    /* getter */
    int getNodeNum() {
        return nodeNum;
    }

    NodeType getNodeType() {
        return nodeType;
    }

    char *getNodeName() {
        return nodeName; 
    }

    Node *getNext() {
        return next;
    }

    Connections *getConList() {
        return conList;
    }

    // PtranFlag getPTranFlag() {
    //     return ptranFlag;
    // }

    /* setter */
    void setNext(Node *nodeIn) {
        next = nodeIn;
    }

    // void setPTranFlag(PtranFlag tmp) {
    //     ptranFlag = tmp;
    // }

    //
    void connect(int conNumIn, Device *deviceIn) {
        consCount++;
        Connections *conPtr = new Connections();
        conPtr->conNum = conNumIn;
        conPtr->device = deviceIn;
        conPtr->next = conList;
        conList = conPtr;
    }

    // 获取连接关系数目
    int getConsCount() {
        return consCount;
    }

    // void printMessage(std::ofstream &outFile);

};


// 节点工厂
class NodeHead {
private:
    NodeHead() {
        nodeList = nullptr;
        nodeCount = 0;
        symtbl.clear();
        id2Name.clear();
        // 默认添加接地节点
        addGroundNode();
    }

    ~NodeHead() {
        Node *node = nodeList, *next;
        while (node != nullptr) {
            next = node->getNext();
            delete node;
            node = next;
        }
        symtbl.clear();
        id2Name.clear();
    }

    // static NodeHead nodeHead;
    Node *nodeList;
    int nodeCount;
    std::unordered_map<std::string, int> symtbl;
    std::unordered_map<int, std::string> id2Name;

public:
    // void addNode(char *nodeName, Device *device, int conNum) {
    //     if (symtbl[nodeName] != 0 || !strcmp(nodeName, "0")) return; // 节点已经存在

    //     // 创建节点
    //     symtbl[nodeName] = ++nodeCount;
    //     Node *tmp_node = new Node(nodeName, nodeCount);
    //     tmp_node->setNext(nodeList);
    //     nodeList = tmp_node;

    //     //维护信息
    //     nodeList->connect(conNum, device);
    //     device->connect(conNum, nodeList);
    // }

    std::string getID2Name(int id) {
        return id2Name[id];
    } 

    void setName2ID(char *name, int id) {
        symtbl[name] = nodeCount + id;
        id2Name[nodeCount + id] = name;
    }

    void addGroundNode();

    void addNode(char *nodeName, Device *device, int conNum, NodeType nodeType);

    int getCount();

    // 获取指定节点编号的节点实例
    Node *getNode(int nodeNum);

    // 获取指定节点名称的节点实例
    Node *getNode(const char *nodeName);

    // 获取头节点
    Node *getNode();

    void printMessage(std::ofstream& outfile);


    // 删除拷贝构造函数和赋值运算符确保单例实例的唯一性
    NodeHead(const NodeHead&) = delete;
    NodeHead& operator=(const NodeHead&) = delete;
    
    static NodeHead& getInstance() {
        static NodeHead nodeHead;
        return nodeHead;
    }

};


/*******************************************NODE  END****************************************************/
/********************************************************************************************************/


/********************************************************************************************************/
/*******************************************Netlist  BEGIN***********************************************/

// 分类类
class Analysis {
private:
    AnalysisType opType, tranType;
    AnalysisType curMode;

    std::vector<std::ostringstream> FOutFile;
    // std::vector<std::ostringstream> JACOutFile;
        
    // // 其它变量
    // std::string outFileName, title;
    // int datum, lastnode, step, total;
    // double tran_initialVal;
    // double tran_step;
    // AnalysisType _type;
    // const int RANDOM_MIN = 333333, RANDOM_MAX = 999999;

    // //
    // double lambda, homostep = 0.05;
    // const double Shomostepratio = 2;
    // const double Fhomostepratio = 0.125;
    // const double minhomostep = 1e-9;
    // const double ERRORGAP = 1e-4;
    // const int ITERATIONNUMS = 20;
    // const int maxSolve = 1e6;
    // std::vector<double> F_X, X, X_n;
    // std::vector<double> G, a;
    // std::vector<std::vector<double>> JAC;

    // tran 分析求解相关
    double tran_stop;
    double tran_step_size;
    const double tran_error_atol = 1e-6;
    const double tran_error_rtol = 1e-6;

    // Ptran 分析求解相关
    double ptran_step_size;
    double pIndctor = 0.00000001;
    double pCapactor = 0.1;
    // const int ptran_max_nriter = 50;
    const int ptran_max_nriter = 10;
    const double ptran_end_time = 1e30;
    const double rtol_ptran = 1e-4;
    const double atol_ptran = 1e-4;
    const double min_ptran_step_size = 1e-9;
    const double max_ptran_step_size = 1e20;
    double init_ptran_step_size = 1.1;

    // NR
    const double rtol_nr = 1e-4;
    const double atol_nr = 1e-4;
    const int maxNRIterNum = 1000;

    //
    int cirNodeTotal;
    int cirBranchVsource;
    int cirBranchInductor;

    const int datum = 0;

public:
    Analysis() {
        opType = AnalysisType::DC_NR;
        tranType = AnalysisType::no_An;
        curMode = opType;

        // ic.clear();
        // nodeset.clear();

        // JACOutFile.resize(cirNodeTotal + cirBranchInductor + cirBranchVsource);

        tran_step_size = ptran_step_size = 1e-9;
    }
    
    ~Analysis() {}

    void setAnalysisType(AnalysisType type) {
        switch (type) {
        case AnalysisType::no_An:
        case AnalysisType::TRAN:
            tranType = type;
            break;
        case AnalysisType::DC_NR:
        case AnalysisType::DC_PTran:
        case AnalysisType::DC_Homotopy:
            opType = type;
            break;
        }
        return;
    }

    void setPtranInductor(double ptranInd) {
        pIndctor = ptranInd;
    }

    void setPtranInitStepSize(double stepsize) {
        init_ptran_step_size = stepsize;
    }

    void setPtranCapactor(double ptranCap) {
        pCapactor = ptranCap;
    }

    void setTranStop(double val) {
        tran_stop = val;
    }

    void setCirNodeTotal(int totalIn) {
        cirNodeTotal = totalIn;
    }

    void setCirBranchVsource(int VSourceNum) {
        cirBranchVsource = VSourceNum;
    }

    void setCirBranchInductor(int InductorNum) {
        cirBranchInductor = InductorNum;
    }

    void setCurMode(AnalysisType typein) {
        curMode = typein;
    }

    void setPtanStepSize(double stepsize) {
        ptran_step_size = stepsize;
    }

    double getInitPtranStepSize() {
        return init_ptran_step_size;
    }

    int getPtranMaxNRIter() {
        return ptran_max_nriter;
    }

    // double getPtranStepSize() {
    //     return ptran_step_size;
    // }

    double getPtranEndTime() {
        return ptran_end_time;
    }

    double getMinPtranStepSize() {
        return min_ptran_step_size;
    }

    double getMaxPtranStepSize() {
        return max_ptran_step_size;
    }

    // double getInitPtranStepSize() {
    //     return init_ptran_step_size;
    // }

    double getRtolNR() {
        return rtol_nr;
    }

    double getAtolNR() {
        return atol_nr;
    }

    AnalysisType getOpType() {
        return opType;
    }

    AnalysisType getTranType() {
        return tranType;
    }

    double getTranStop() {
        return tran_stop;
    }

    double getPTranStepSize() {
        return ptran_step_size;
    }

    double getPInductor() {
        return pIndctor;
    }

    double getPCapactor() {
        return pCapactor;
    }

    double getPtranAtol() {
        return atol_ptran;
    }

    double getPtranRtol() {
        return rtol_ptran;
    }

    double getTranStepSize() {
        return tran_step_size;
    }

    double getMinTranStepSize() {
        return min_ptran_step_size;
    }

    int getCirNodeTotal() {
        return cirNodeTotal;
    }

    int getCirBranchVsource() {
        return cirBranchVsource;
    }

    int getCirBranchInductor() {
        return cirBranchInductor;
    }

    int getDatum() {
        return datum;
    }

    int getMaxNRIterNum() {
        return maxNRIterNum;
    }

    AnalysisType getCurMode() {
        return curMode;
    }

    std::vector<std::ostringstream>& getFOutFileAndClear() {
        FOutFile.clear();
        FOutFile.resize(cirNodeTotal + cirBranchVsource + cirBranchInductor);
        return FOutFile;
    }

};

// 网表类
class Netlist { 
private:
    ModelHead &modelList = ModelHead::getInstance();
    DeviceHead &deviceList = DeviceHead::getInstance();
    NodeHead &nodeList = NodeHead::getInstance();
    // Analysis *analysis;

    char title[NameLength];
    const int datum = 0;

    int COOopNum;
    SMPMatrix *matrix;
    COOElement *linkedCOOList;
    std::vector<std::ostringstream> FOutFile;

    double tran_stop;

    std::unordered_map<std::string, double> ic;
    std::unordered_map<std::string, double> nodeset;

public:
    Netlist() {
        // analysis = new Analysis();
        matrix = nullptr;
        linkedCOOList = nullptr;
        COOopNum = 0;
    }

    ~Netlist() {
        // delete analysis;
        // destroySMPMatrix(matrix);
    }

    void loadIC(double *X) {
        for (auto u:ic) {
            int id = nodeList.getNode(u.first.c_str())->getNodeNum() - 1;
            X[id] = u.second;
        }
    }

    void loadNodeSet(double *X) {
        for (auto u:nodeset) {
            Node *node = nodeList.getNode(u.first.c_str());
            if (node == nullptr) continue;
            int id = node->getNodeNum() - 1;
            X[id] = u.second;
        }
    }

    void setIC(std::string nodeName, double value) {
        ic[nodeName] = value;
    }

    void setNodeSet(std::string nodeName, double value) {
        nodeset[nodeName] = value;
    }

    void setTitle(char *nameIn) {
        assert(strlen(nameIn) <= NameLength);
        strcpy(title, nameIn);
    }

    // Analysis* getAnalysis() {
    //     return analysis;
    // }
    
    ModelHead& getModelHead() {
        return modelList;
    }

    DeviceHead& getDeviceHead() {
        return deviceList;
    }

    NodeHead& getNodeHead() {
        return nodeList;
    }

    char* getTitle() {
        return title;
    }

    int getDatum() {
        return datum;
    }

};

/*******************************************Netlist  END*************************************************/
/********************************************************************************************************/

#endif
