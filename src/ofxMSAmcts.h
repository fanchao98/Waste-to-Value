/*
一个非常简单的基于C++11的模板化MCTS（蒙特卡洛树搜索）实现，附带openFrameworks示例。

MCTS代码基于Java（Simon Lucas - University of Essex）和Python（Peter Cowling，Ed Powley，Daniel Whitehouse - University of York）的实现，网址：http://mcts.ai/code/index.html
*/

#pragma once

#include "TreeNodeT.h"
#include "MSALoopTimer.h"
#include <cfloat>
#include"IState.h"
#include"ReFab.h"

namespace msa {
    namespace mcts {

        // 状态必须符合状态接口（参见IState.h）
        // 动作可以是任何事物（你的State类知道如何处理）
        //template <class State, typename Action>
        class UCT {
            typedef TreeNodeT TreeNode;

        private:
            LoopTimer timer;
            int iterations;
            
        public:
            ReFab refab;
            float uct_k;					// UCT函数中的k值。默认值= sqrt(2)
            unsigned int max_iterations;	// 最多执行这么多次迭代（0表示一直运行）
            unsigned int max_millis;		// 最多运行这么多毫秒（0表示一直运行）
            unsigned int simulation_depth;	// 运行模拟的帧数

            //--------------------------------------------------------------
            UCT() :
                iterations(0),
                uct_k(sqrt(2)),  //sqrt(2)  //0.4
                max_iterations(500),  //1000  //100 //1500
                max_millis(0),
                simulation_depth(50)
            {}

            /*Get_ReFab(ReFab& refab) {
				this->refab = refab;
            }*/

            //--------------------------------------------------------------
            const LoopTimer& get_timer() const {
                return timer;
            }

            const int get_iterations() const {
                return iterations;
            }


            //--------------------------------------------------------------
            // 根据uct分数获取给定TreeNode的最佳（直接的）子节点
            TreeNode* get_best_uct_child(TreeNode* node, float uct_k) const {
                // 检查合法性
                if (!node->is_fully_expanded()) return NULL;

                float best_utc_score = -std::numeric_limits<float>::max();
                TreeNode* best_node = NULL;

                // 遍历所有直接子节点并找到最佳的UTC分数，若是叶子节点，则换一个
                int num_children = node->get_num_children();
                for (int i = 0; i < num_children; i++) {
                    TreeNode* child = node->get_child(i);
                    float uct_exploitation = (float)child->get_value() / (child->get_num_visits() + FLT_EPSILON);
                    float uct_exploration = sqrt(log((float)node->get_num_visits() + 1) / (child->get_num_visits() + FLT_EPSILON));
                    float uct_score = uct_exploitation + uct_k * uct_exploration;


                    if (uct_score > best_utc_score) {
                        best_utc_score = uct_score;
                        best_node = child;
                    }
                }

                return best_node;
            }


            //--------------------------------------------------------------
            TreeNode* get_most_visited_child(TreeNode* node) const {
                int most_visits = -1;
                TreeNode* best_node = NULL;

                // 遍历所有直接子节点并找到访问次数最多的
                int num_children = node->get_num_children();
                for (int i = 0; i < num_children; i++) {
                    TreeNode* child = node->get_child(i);
                    if (child->get_num_visits() > most_visits) {
                        most_visits = child->get_num_visits();
                        best_node = child;
                    }
                }

                return best_node;
            }



            //--------------------------------------------------------------
            Action run(State& current_state, int interleaving_iteration, unsigned int seed = 1, vector<State>* explored_states = nullptr) {
                // 初始化计时器
                timer.init();


                // 使用当前状态初始化根TreeNode
                int cont_node_id = 0;
                vector<vector<Point_3>> all_vertices_in_components;
                current_state.total_num_of_planes = refab.Calculate_current_num_of_connected_components(current_state.V_offset, current_state.V_ori, current_state.total_num_of_complete_planes,
                    current_state.barycenter_of_planes,current_state.selected_planes_normal, all_vertices_in_components,interleaving_iteration);
                current_state.expected_normals_of_cutting_planes = refab.expected_normals_of_cutting_planes;
                current_state.mean_scalar_value_of_cutting_planes = refab.mean_scalar_value_of_cutting_planes;
                current_state.additive_interface_points = refab.additive_interface_points;
                current_state.index_of_additive_interface = refab.index_of_additive_interface;
                TreeNode root_node(current_state);
                root_node.node_id = cont_node_id;
                TreeNode* best_node = NULL;
                
                // 迭代
                iterations = 1;
                double sum_volume = 0;
                int cont_right = 0;
                double max_modified_volume = 0;
                double last_max_volume = 0;
                double new_volume = 0;
                ofstream file_save_volume("MCTS_temp\\" + refab.mesh_target + "\\save_volume.txt");
                while (true) {
                    // 指示循环开始
                    timer.loop_start();

                    clock_t start_time_total, end_time_total;
                    start_time_total = clock();
                    // 1. 选择。从根开始，在所有完全展开的节点上使用UCT挖掘树
                    TreeNode* node = &root_node;
                    while (!node->is_terminal_2() && node->is_fully_expanded()) {
                        node = get_best_uct_child(node, uct_k);
                        //	assert(node);	// 检查合法性
                    }

                    // 2. 扩展，通过添加单个子节点（如果不是终端或没有完全展开）
                    if (!node->is_fully_expanded() && !node->is_terminal()) {
                        if(node->get_parent()!=NULL)
                            node->new_current_points = node->get_parent()->new_current_points;   //从父节点获得更新后的移除点集
                        node = node->expand(all_vertices_in_components, refab.cont_components_of_each_planes[0]);
                        cont_node_id++;
                        node->node_id = cont_node_id;
                    }
                        

                    State state(node->get_state());          
                    state.expected_normals_of_cutting_planes = refab.expected_normals_of_cutting_planes;
                    state.mean_scalar_value_of_cutting_planes = refab.mean_scalar_value_of_cutting_planes;
                    state.additive_interface_points = refab.additive_interface_points;
                    state.index_of_additive_interface = refab.index_of_additive_interface;

                    // 3. 模拟
                    //if (state.is_terminal()) break;
                    bool flag_satisfy_constraints;
                    vector<Plane> parent_planes;
                    parent_planes.clear();
                    for (int i = 0; i <state.normal_of_planes.size(); i++)
                    {
						Plane plane(state.normal_of_planes[i], state.origin_of_planes[i]);
						parent_planes.push_back(plane);
					}

                    clock_t start_time, end_time;
                    start_time = clock();
                    vector<Plane> best_modifed_cutting_planes(0);
                    vector<Point_3> new_current_removed_points;
                    bool is_max_volume;
                    refab.saved_removed_triangles.load("MCTS_temp\\" + refab.mesh_target + "\\removed_triangles-" + to_string(0) + ".obj");
                    refab.saved_new_D_I.load("MCTS_temp\\" + refab.mesh_target + "\\new_D_I-" + to_string(0) + ".obj");
                    flag_satisfy_constraints = refab.Modify_Cutting_Planes_2(parent_planes, new_volume, best_modifed_cutting_planes,max_modified_volume, last_max_volume, new_current_removed_points,false, is_max_volume);
                    end_time = clock();
                    std::cout << "Simulation time: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;
                    if (flag_satisfy_constraints == true)
                        cont_right++;
                    sum_volume += new_volume;

                    // 获取所有代理的奖励向量
                    const std::vector<float> rewards = state.evaluate(flag_satisfy_constraints, new_volume, is_max_volume, last_max_volume);
                    //cout <<  30 * (state.V_current - state.V_offset) / (state.V_ori - state.V_offset);  

                    // 添加到历史记录
                    if (explored_states) explored_states->push_back(state);

                    node->new_current_points = new_current_removed_points;
                    // 4. 反向传播
                    cout << "current node id: " << node->node_id << endl;
                    while (node) {
                        node->update(rewards);
                        node = node->get_parent();
                        if(node != NULL)
                            cout <<"parent node id: " << node->node_id << endl;
                    }

                    // 找到访问最多的子节点
                    best_node = get_most_visited_child(&root_node);

                    // 指示计时器的循环结束
                    timer.loop_end();

                    // 如果当前总运行时间（自初始化以来）超过max_millis，则退出循环
                    if (max_millis > 0 && timer.check_duration(max_millis)) break;

                    // 如果当前迭代次数超过max_iterations，则退出循环
                    if (max_iterations > 0 && iterations > max_iterations) break;

                    cout << "Iteration: " << iterations << endl;
                    if (rewards[0] != 0)
                        cout << "A √ nodes!  Reward: " << rewards[0] <<endl << endl;
                    else
                        cout << "A × nodes." << endl << endl;
                    file_save_volume << new_volume/10000.0 <<endl;
                    iterations++;

                    end_time_total = clock();
                    //std::cout << "once iteration time: " << (double)(end_time_total - start_time_total) / CLOCKS_PER_SEC << "s" << std::endl;
                }

                cout << "total √ nodes: " << cont_right << endl;
                cout << "sum volume: " << sum_volume << endl;
                vector<Plane> best_modifed_cutting_planes = refab.modified_cutting_planes_by_MTCS;
                refab.record_data[6] = best_modifed_cutting_planes.size();
                vector<Plane> parent_planes;
                vector<Point_3> new_current_removed_points;
                bool is_max_volume;
                refab.saved_removed_triangles.load("MCTS_temp\\" + refab.mesh_target + "\\removed_triangles-" + to_string(0) + ".obj");
                refab.saved_new_D_I.load("MCTS_temp\\" + refab.mesh_target + "\\new_D_I-" + to_string(0) + ".obj");
                refab.Modify_Cutting_Planes_2(parent_planes, new_volume, best_modifed_cutting_planes, max_modified_volume, last_max_volume,new_current_removed_points,true, is_max_volume);
                cout << "Final volume: " << new_volume << endl;
                refab.record_data[7] = new_volume / refab.volume_of_target_model;

                //更新 Selected_Planes，用于下一次迭代
                refab.Selected_Planes = best_modifed_cutting_planes;
                refab.flag_useful_cutting_planes = new bool[best_modifed_cutting_planes.size()];
                for(int t =0;t < best_modifed_cutting_planes.size();t++)
					refab.flag_useful_cutting_planes[t] = true;

                // 返回最佳节点的动作
                if (best_node) return best_node->get_action();

                // 我们不应该到达这里
                return Action();
            }


        };
    }
}
