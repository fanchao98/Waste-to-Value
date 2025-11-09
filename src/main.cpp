#include "ReFab.h"
#include "IState.h"
#include "ofxMSAmcts.h"
#include "PPS.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <filesystem>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <direct.h>  // for _mkdir

using namespace msa::mcts;
using namespace std;

int main(int argc, char* argv[])
{
    clock_t start_time1, start_time2, start_time3, start_time4, start_time5;
    clock_t end_time1, end_time2, end_time3, end_time4, end_time5;

    // Tool and nozzle parameters
    nozzle Nozzle;
    Nozzle.upper_surface_r = 4.5;   // upper surface radius
    Nozzle.lowwer_surface_r = 2.0;  // lower surface radius
    Nozzle.nozzle__H_total = 10;    // total nozzle height
    Nozzle.nozzle_H_half = 5;       // half nozzle height

    cutter Cutter;
    Cutter.cylinder_height = 13;    // cutter height
    Cutter.ball_r = 1.5;            // ball radius
    Cutter.carriage_r = 25;         // carriage radius
    Cutter.carriage_height = 30;    // carriage height

    // Number of interleaving iterations (currently unused)
    int threshold_interleaving_iteration = 1;

    // Set root path (relative)
    std::string rootpath = "./";
    std::string mesh_ori = "F_coral_ori";
    std::string mesh_target = "F_coral_target";

    // Create ReFab object
    ReFab refab(rootpath, Cutter, Nozzle);
    refab.record_data.resize(8);

    refab.open_vis_green_points = false;
    refab.open_vis_red_points = true;

    // Input meshes
    std::string meshname1 = "models/" + mesh_ori + ".off";
    std::string meshname2 = "models/" + mesh_target + ".off";
    std::string intersection = "output/" + mesh_target + "/D_I.off";
    refab.mesh_target = mesh_target;

    // Create necessary folders
    string path = refab.rootPath + "output/" + mesh_target;
    _mkdir(path.c_str());
    path = refab.rootPath + "temp_vis/" + mesh_target;
    _mkdir(path.c_str());
    path = refab.rootPath + "MCTS_temp/" + mesh_target;
    _mkdir(path.c_str());

    ofstream out_record(refab.rootPath + "output/" + mesh_target + "/record.txt");

    start_time1 = clock();
    //refab.ReOrientation(meshname1, meshname2, 0);
    end_time1 = clock();
    meshname2 = meshname2.substr(0, meshname2.size() - 4) + "_new.off";
    out_record << "Orientation time: " << (double)(end_time1 - start_time1) / CLOCKS_PER_SEC << "s" << std::endl;

    // ------------------- Mesh Boolean Operations -------------------
    start_time2 = clock();
    refab.Intersection_and_Split(meshname1, meshname2,
        "output/" + mesh_target + "/Interface_S.off",
        "output/" + mesh_target + "/Interface_A.off",
        "output/" + mesh_target + "/D_I.off", 0);

    refab.Difference_two_mesh(meshname1, meshname2, "output/" + mesh_target + "/D_S.off", 0);
    refab.Difference_two_mesh(meshname2, meshname1, "output/" + mesh_target + "/D_A.off", 0);
    refab.Modify_DS_and_DA("output/" + mesh_target + "/D_S.off",
        "output/" + mesh_target + "/D_A.off",
        "output/" + mesh_target + "/hang_in_area.off");
    refab.DetermineBase(meshname1);

    // ------------------- Call Python Poisson Sampling Script -------------------
    std::string command = "python poisson_sample.py " + mesh_target;
    int ret = std::system(command.c_str());
    if (ret == 0) {
        std::cout << "Python script executed successfully!\n";
    }
    else {
        std::cerr << "Python script failed with code: " << ret << std::endl;
    }

    // ------------------- Sampling and Collision Detection -------------------
    refab.Mesh_Regular_Sampling_C1(intersection, 0.03);
    refab.Get_sampling_points_in_DS_and_DA("output/" + mesh_target + "/D_S.xyzn",
        "output/" + mesh_target + "/D_A.xyzn");
    refab.Collision_Detection_For_DS_and_DA();
    end_time2 = clock();

    // ------------------- Initialize Cutting Planes -------------------
    start_time3 = clock();
    refab.Initialize_Cutting_Planes("output/" + mesh_target + "/Interface_A.off");
    end_time3 = clock();

    // Record statistics
    out_record << "RV: " << refab.record_data[0] << endl;
    out_record << "A: " << refab.record_data[1] << endl;
    out_record << "PS: " << refab.record_data[2] << endl;
    out_record << "PA: " << refab.record_data[3] << endl;

    // ------------------- Run Monte Carlo Tree Search -------------------
    start_time4 = clock();
    int interleaving_iteration = 0;
    UCT* ucts = new UCT[threshold_interleaving_iteration];
    State state;
    state.reset();
    UCT uct = ucts[interleaving_iteration];
    uct.refab = refab;
    uct.run(state, interleaving_iteration); // run MCTS
    end_time4 = clock();

    // ------------------- Output Results -------------------
    out_record << "*----------------------------------------------------*" << endl;
    out_record << "Pre-process time: " << (double)(end_time2 - start_time2) / CLOCKS_PER_SEC << "s" << std::endl;
    out_record << "Initialization time: " << (double)(end_time3 - start_time3) / CLOCKS_PER_SEC << "s" << std::endl;
    out_record << "MCTS time: " << (double)(end_time4 - start_time4) / CLOCKS_PER_SEC << "s" << std::endl;
    out_record << "*----------------------------------------------------*" << endl;

    // Record performance metrics
    out_record << "RV: " << refab.record_data[0] << endl;
    out_record << "A: " << refab.record_data[1] << endl;
    out_record << "PS: " << refab.record_data[2] << endl;
    out_record << "PA: " << refab.record_data[3] << endl;
    out_record << "RI: " << uct.refab.record_data[4] << endl;
    out_record << "Cb: " << uct.refab.record_data[5] << endl;
    out_record << "Ca: " << uct.refab.record_data[6] << endl;
    out_record << "RF: " << uct.refab.record_data[7] << endl;

    return 0;
}
