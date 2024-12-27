#include "minispice.hpp"

#include <vector>
#include <string>
#include <Python.h>

int plot(std::unique_ptr<Analysis>& analysis, std::unique_ptr<Netlist>& netlist) {
    
    std::vector<std::vector<double>>& results = analysis->getPtranResults();
    std::vector<double>& timepoints = analysis->getPtranTimepoints();
    std::unordered_map<int, std::string>& id2Name = netlist->getNodeHead().getID2Name_();

    // 
    Py_Initialize();

    // 屏蔽警告
    PyRun_SimpleString("import warnings\n"
                       "warnings.filterwarnings('ignore', category=UserWarning)");

    if (!Py_IsInitialized()) {
        std::cerr << "Python initialization failed!" << std::endl;
        return MINISPICEERROR;
    }

    // 添加路径
    std::string script_path = PLOT_SCRIPT_PATH;
    PyObject* sys_path = PySys_GetObject("path"); // 获取 sys.path
    PyObject* py_script_dir = PyUnicode_FromStringAndSize(script_path.c_str(), script_path.size());
    PyList_Append(sys_path, py_script_dir); // 添加脚本目录到 sys.path
    Py_DECREF(py_script_dir);

    // 引入绘图脚本
    PyObject* pModule = PyImport_ImportModule("plot");
    if (!pModule) {
        std::cerr << "Failed to load Python module 'plot'" << std::endl;
        Py_Finalize();
        return MINISPICEERROR;
    }

    // 获取Python函数
    PyObject* pFunc = PyObject_GetAttrString(pModule, "visualize");
    if (!pFunc || !PyCallable_Check(pFunc)) {
        std::cerr << "Python function 'visualize' not found or not callable" << std::endl;
        Py_DECREF(pModule);
        Py_Finalize();
        return MINISPICEERROR;
    }

    // 将 std::vector 转换为 Python 的 list
    PyObject* pListResults = PyList_New(results.size());
    for (int i = 0; i < results.size(); ++i) {
        PyObject* pInnerList = PyList_New(results[i].size());
        for (int j = 0; j < results[i].size(); ++j) {
            PyList_SetItem(pInnerList, j, PyFloat_FromDouble(results[i][j]));  // 设置节点电压变化曲线值
        }
        PyList_SetItem(pListResults, i, pInnerList);  // 将内层列表加入外层列表
    }

    PyObject* pListTimepoints = PyList_New(timepoints.size());
    for (int i = 0; i < timepoints.size(); ++i) {
        PyList_SetItem(pListTimepoints, i, PyFloat_FromDouble(timepoints[i]));  // 设置节点电压变化时间值
    }

    PyObject* pListNodeName = PyList_New(results[0].size());
    for (int i = 0; i < results[0].size(); ++i) {
        PyList_SetItem(pListNodeName, i, PyUnicode_FromStringAndSize(id2Name[i + 1].c_str(), id2Name[i + 1].size()));  // 设置节点名
    }
    
    // 调用函数
    PyObject* pArgs = PyTuple_New(3);
    PyTuple_SetItem(pArgs, 0, pListResults);
    PyTuple_SetItem(pArgs, 1, pListTimepoints);
    PyTuple_SetItem(pArgs, 2, pListNodeName);

    PyObject* pRetValue = PyObject_CallObject(pFunc, pArgs);

    //
    if (!pRetValue) {
        PyErr_Print();  // 打印 Python 错误信息
        std::cerr << "Error executing Python function." << std::endl;
    }

    // 清空PyObject
    Py_DECREF(pListResults);
    Py_DECREF(pListTimepoints);
    Py_DECREF(pListNodeName);
    Py_DECREF(pArgs);
    Py_DECREF(pRetValue);
    Py_DECREF(pFunc);
    Py_DECREF(pModule);

    // 终止Py环境
    Py_Finalize();

    return MINISPICEOK;
}