#include "minispice.hpp"

#include <vector>
#include <string>
// #include <Python.h>

int plot(std::unique_ptr<Analysis>& analysis) {
    
//     // 很重要
//     Py_SetPythonHome(L"D:/Anaconda3");
    
//     // 初始化python环境
//     Py_Initialize();
//     if (!Py_IsInitialized()) {
//         printf("Py_Initialize failed!!!\n");
//         return -1;
//     }

//     // 测试语句
//     // PyRun_SimpleString("print('Hello Python!')\n");
//     PyRun_SimpleString("import os,sys");//执行import语句，把当前路径的上一级加入路径中，为了找到plot.py
//     PyRun_SimpleString("sys.path.append('../scripts')");
//     // PyRun_SimpleString("print(os.getcwd())");//测试打印当前路径

//     // 定义调用函数时的相关变量
//     PyObject *pModule;
//     PyObject *pFunction;
//     PyObject *pArgs;
//     PyObject *pRetValue;

//     // 加载脚本进来
//     pModule = PyImport_ImportModule("plot"); // 注意这里一定要是plot，不带.py
//     if (!pModule) {
//         printf("import python failed!!!\n");
//         return -1;
//     }

//     // 查找工具函数
//     pFunction = PyObject_GetAttrString(pModule, "plot_util");
//     if (!pFunction) {
//         printf("get python function failed!!!\n");
//         return -1;
//     }

//     // 构建参数
//     // 预处理vector(cpp) -> list(py)
//     PyObject *list1 = PyList_New(X.size()); // X -> x_axis_data
//     for (int i = 0; i < X.size(); i++) {
//         PyList_SetItem(list1, i, Py_BuildValue("d", X[i]));
//     }
//     PyObject *list2 = PyList_New(0);
//     for (int i = 0; i < Y.size(); i++) {
//         PyObject *t = PyList_New(Y[i].size());
//         for (int j = 0; j < Y[i].size(); j++) {
//             PyList_SetItem(t, j, Py_BuildValue("d", Y[i][j]));
//         }
//         PyList_Append(list2, t);
//         Py_DECREF(t);   
//     }
    
//     pArgs = PyTuple_New(6);
//     PyTuple_SetItem(pArgs, 0, list1);
//     PyTuple_SetItem(pArgs, 1, list2);
//     PyTuple_SetItem(pArgs, 2, Py_BuildValue("s", x_name.c_str()));
//     PyTuple_SetItem(pArgs, 3, Py_BuildValue("s", y_name.c_str()));
//     PyTuple_SetItem(pArgs, 4, Py_BuildValue("s", title.c_str()));
//     PyTuple_SetItem(pArgs, 5, Py_BuildValue("s", ("fig" + std::to_string(step)).c_str()));

//     // 调用函数
//     pRetValue = PyObject_CallObject(pFunction, pArgs);

//     // 清空PyObject
//     Py_DECREF(pModule);
//     Py_DECREF(pFunction);
//     Py_DECREF(pArgs);
//     Py_DECREF(pRetValue);
//     Py_DECREF(list1);
//     Py_DECREF(list2);

//     // 终止Py环境
//     Py_Finalize();

//     return 512;

    return MINISPICEOK;
}