#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include "minispice.hpp"
#include "device.hpp"
#include "config.hpp"

int main(int argc, char** argv) {

    // 命令行参数处理
    parser_cmd::parse_args(argc, argv);

    // 
    MPI_Init(&argc, &argv);
    int world_size, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#if !defined(MINISPICESUPERLU) && !defined(MINISPICEMUMPS)
    if (world_size > 1) {
        if (myid == 0) {
            std::cout << "KLU, Eigen and MiniSPICELU don't need multiprocess env!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
#endif

    //
    if (myid == 0) {
        if(parser_cmd::instrError) {
            std::cout << argv[0] << " Invalid instruction, please use '-h' to understand the usage!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        } else if(parser_cmd::help) {
            std::cout << "Usage: " << argv[0] << "<-h> -f FILENAME <-o OUTPUT_FILE>" << std::endl;
            std::cout << "  ------- Listing options -------" << std::endl;
            std::cout << "  -h              Print usage and this help message." << std::endl;
            std::cout << "  -o[- filename]  When is '-', it means the output stream is standard output; When 'filename', represents the output stream as the file." << std::endl;
            std::cout << "  -f[ filename ]  'filename' represents the input stream as the file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //
        if (parser_cmd::inFileName.empty()) {
            std::cout << "Please input netlist file!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // 主要变量初始化
    std::unique_ptr<Netlist> netlist = std::make_unique<Netlist>();
    std::unique_ptr<Analysis> analysis = std::make_unique<Analysis>();
    std::unique_ptr<SMPMatrix> matrix = std::make_unique<SMPMatrix>();

#ifdef MINISPICESUPERLU
    // 确定nprow和npcol的大小【列优先大】
    int nprow, npcol;
    switch (world_size) {
        case 1:
            npcol = 1;
            nprow = 1;
            break;
        case 4:
            npcol = 2;
            nprow = 2;
            break;
        case 8:
            npcol = 4;
            nprow = 2;
            break;
        case 16:
            npcol = 4;
            nprow = 4;
            break;
        case 32:
            npcol = 8;
            nprow = 4;
            break;
        case 64:
            npcol = 8;
            nprow = 8;
            break;
        default:
            if (myid == 0) {
                std::cout << "Invalid processor number, please <1/2/4/8/16/32/64> usage!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
    }
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &matrix->superluMatrix->grid);
    // std::cout << myid << " s " <<  matrix->superluMatrix->grid.comm << " s " << std::endl;
#endif

#ifdef MINISPICESUPERLU
    myid = matrix->superluMatrix->grid.iam;
#else
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif

#ifdef MINISPICEMUMPS
    matrix->mumpsMatrix->myid = myid;
#endif

MPI_ENV_0_PID_BEGIN(myid)
    
    std::cout << "=====================Simulation Start=====================" << std::endl;

    std::cout << "Stage1: Parsing Netlist Start" << std::endl;

    // 1. 解析网表文件，存储信息到nodeList、compList、modelList中
    int error = parser(netlist, analysis, parser_cmd::inFileName, parser_cmd::outFileName);

    if (error == MINISPICEERROR) {
        std::cerr << "Stage1: Parser netlist failed!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    } else {
        std::cout << "Stage1: Parser netlist completed!" << std::endl;
    }
    
    // // 2. setup
    // error = setup(matrix, netlist, analysis, parser_cmd::outFileName);
    
    // if (error == MINISPICEERROR) {
    //     std::cerr << "Stage2: Setup failed!" << std::endl;
    //     MPI_Abort(MPI_COMM_WORLD, 1);
    // }

    std::cout << "Stage2: Dojob Start" << std::endl;

MPI_ENV_0_PID_END

    // 3. dojob
    int error = dojob(matrix, netlist, analysis, parser_cmd::outFileName, myid);

MPI_ENV_0_PID_BEGIN(myid)

    if (error == MINISPICEERROR) {
        std::cerr << "Stage2: Dojob failed!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    } else {
        std::cout << "Stage2: Dojob completed!" << std::endl;
    }

    // 4. plot
    error = plot(analysis);

    if (error == MINISPICEERROR) {
        std::cerr << "Stage4: Visualization failed!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::cout << "=====================Simulation Final=====================" << std::endl;

MPI_ENV_0_PID_END


#ifdef MINISPICESUPERLU
    superlu_gridexit(&matrix->superluMatrix->grid);
#endif

    MPI_Finalize();

    return 0;
}