#include <sys/time.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <utility>
#include <vector>
#include <iostream>
#include <filesystem>
#include <string>
#include <queue>      
#include <climits>
#include <unordered_set>
#include <cassert>
#include <stack>
#include <numeric>
#include <unordered_map>
#include <sstream>
#include <map>
#include <functional>
#include <set>
#include <cmath>
#include <random>


namespace bs {

using namespace std;

struct node {
  int id;              
  int N_O_SZ, N_I_SZ;
  int *N_O, *N_I;
  int vis;  
	
  int group_id;        // 该节点所属的组ID，-1表示未分组（汇点）
	
	
    
 vector<vector<pair<int, int>>> interval_array;
   int tin, tout;                        
   int t;                         
};

vector<node> nodes;
vector<vector<int>> node_groups;  // node_groups[i] 存储第i组的所有节点ID

int vis_cur, cur;
int subgraph_number;
vector<int> group_sizes;

static constexpr int BS_THRESHOLD = 8;


void read_graph(const char *filename) {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);

    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "无法打开文件 %s\n", filename);
        return;
    }

    char buffer[1024];
    // 跳过第一行注释
    if (!fgets(buffer, sizeof(buffer), file)) {
        fprintf(stderr, "无法读取注释行\n");
        fclose(file);
        return;
    }

    int n, m;
    // 读第二行：节点数 n 和 边数 m
    if (fscanf(file, "%d %d", &n, &m) != 2) {
        fprintf(stderr, "读取 n, m 失败\n");
        fclose(file);
        return;
    }

    nodes.resize(n);
    vector<vector<int>> N_O(n), N_I(n);

    for (int i = 0; i < n; i++) {
        int u;
        if (fscanf(file, "%d", &u) != 1) {
            break;
        }

        int c;
        while ((c = fgetc(file)) != EOF && c != '\n') {
            if (isdigit(c) || c == '-') {
                ungetc(c, file);
                int v;
                if (fscanf(file, "%d", &v) == 1) {
                    N_O[u].push_back(v);
                    N_I[v].push_back(u);
                }
            }
        }
    }

    fclose(file);

    for (int u = 0; u < n; u++) {
        nodes[u].id = u;         // 设置节点ID
				
        nodes[u].N_O_SZ = (int)N_O[u].size();
        nodes[u].N_O = new int[nodes[u].N_O_SZ];
        for (int i = 0; i < nodes[u].N_O_SZ; i++)
            nodes[u].N_O[i] = N_O[u][i];

        nodes[u].N_I_SZ = (int)N_I[u].size();
        nodes[u].N_I = new int[nodes[u].N_I_SZ];
        for (int i = 0; i < nodes[u].N_I_SZ; i++)
            nodes[u].N_I[i] = N_I[u][i];

        // 初始化 interval_array，大小为 subgraph_number
        nodes[u].interval_array.resize(subgraph_number);
        
        // 初始化其他字段
        nodes[u].vis = 0;
        nodes[u].group_id = -1;  // -1 表示未分组
        nodes[u].tin = -1;
        nodes[u].tout = -1;
    }

    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;
    printf("read time(graph): %.3f ms\n", elapsed * 1000);
}


void calculate_subgraph_sizes() {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);
    
    printf("\n=== 计算子图分组大小（随机平均分配所有节点） ===\n");
    
    int total_nodes = nodes.size();
    printf("总节点数: %d\n", total_nodes);
    printf("需要分成 %d 个组\n", subgraph_number);
    
    if (subgraph_number <= 0) {
        printf("错误：subgraph_number必须大于0\n");
        return;
    }
    
    if (total_nodes == 0) {
        printf("错误：没有节点可以分组\n");
        return;
    }
    
    // 1. 初始化数据结构
    group_sizes.resize(subgraph_number);
    node_groups.resize(subgraph_number);
    
    // 2. 计算每组的大小
    int base_size = total_nodes / subgraph_number;
    int remainder = total_nodes % subgraph_number;
    
    for (int i = 0; i < subgraph_number; i++) {
        if (i < remainder) {
            group_sizes[i] = base_size + 1;
        } else {
            group_sizes[i] = base_size;
        }
        node_groups[i].reserve(group_sizes[i]);  // 预分配空间
    }
    
    // 3. 创建随机排列的节点ID列表
    vector<int> node_ids(total_nodes);
    for (int i = 0; i < total_nodes; i++) {
        node_ids[i] = i;
    }
    
    
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::shuffle(node_ids.begin(), node_ids.end(), gen);
    
    // 4. 将节点分配到各组
    int node_index = 0;
    for (int group = 0; group < subgraph_number; group++) {
        for (int i = 0; i < group_sizes[group]; i++) {
            int node_id = node_ids[node_index++];
            nodes[node_id].group_id = group;
            node_groups[group].push_back(node_id);
        }
    }
    
    // 5. 输出分组结果
    printf("\n分组结果（随机平均分配）：\n");
    int total_check = 0;
    int min_size = INT_MAX;
    int max_size = 0;
    
    for (int i = 0; i < subgraph_number; i++) {
        printf("  组 %d: %d 个节点\n", i, group_sizes[i]);
        total_check += group_sizes[i];
        min_size = min(min_size, group_sizes[i]);
        max_size = max(max_size, group_sizes[i]);
    }
    
    // 6. 验证
    printf("\n验证：\n");
    printf("  总分配节点数: %d\n", total_check);
    printf("  应分配节点数: %d\n", total_nodes);
    printf("  最小组大小: %d\n", min_size);
    printf("  最大组大小: %d\n", max_size);
    printf("  大小差异: %d (应该 <= 1)\n", max_size - min_size);
    
    
    
    // 8. 验证所有节点都被分配
    int unassigned_count = 0;
    for (int i = 0; i < total_nodes; i++) {
        if (nodes[i].group_id == -1) {
            unassigned_count++;
        }
    }
    if (unassigned_count > 0) {
        printf("  警告：有 %d 个节点未被分配！\n", unassigned_count);
    }
    
    // 9. 计算标准差（衡量分配的均匀程度）
    double mean = (double)total_nodes / subgraph_number;
    double variance = 0;
    for (int i = 0; i < subgraph_number; i++) {
        variance += pow(group_sizes[i] - mean, 2);
    }
    variance /= subgraph_number;
    double std_dev = sqrt(variance);
    
    printf("  平均每组大小: %.2f\n", mean);
    printf("  标准差: %.2f (越小越均匀)\n", std_dev);
    
    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;
    printf("\n计算时间: %.3f ms\n", elapsed * 1000);
}



// 优化的合并区间列表 - 确保结果始终排序
void merge_intervals(vector<pair<int,int>> &to, const vector<pair<int,int>> &from) {
    if (from.empty()) return;
    if (to.empty()) {
        to = from;
        return;
    }
    
    // 预先reserve空间
    to.reserve(to.size() + from.size());
    to.insert(to.end(), from.begin(), from.end());
    
    // 排序所有区间
    sort(to.begin(), to.end());
    
    // 原地合并重叠区间
    int write_idx = 0;
    for (int read_idx = 1; read_idx < to.size(); ++read_idx) {
        if (to[write_idx].second >= to[read_idx].second) {
            // 当前区间被包含，跳过
            continue;
        } else if (to[write_idx].second >= to[read_idx].first - 1) {
            // 重叠或相邻，合并
            to[write_idx].second = max(to[write_idx].second, to[read_idx].second);
        } else {
            // 不相邻，添加新区间
            to[++write_idx] = to[read_idx];
        }
    }
    to.resize(write_idx + 1);
    
}








void build_ferrari_index_for_subgraph(int group) {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);
    
    printf("构建组 %d 的Ferrari索引\n", group);
    
    int n = nodes.size();
    
    for (int u = 0; u < n; u++) {
        
        nodes[u].tin = -1;
        nodes[u].tout = -1;
    
    }
    // 1. 拓扑排序
      vector<int> indeg(n, 0), topo_order;
      for (int i = 0; i < n; ++i) {
          for (int j = 0; j < nodes[i].N_O_SZ; ++j) {
              indeg[nodes[i].N_O[j]]++;
          }
      }
      
      queue<int> q;
      for (int i = 0; i < n; ++i) {
          if (indeg[i] == 0) {
              q.push(i);
          }
      }

      while (!q.empty()) {
          int v = q.front(); 
          q.pop();
          topo_order.push_back(v);
          for (int i = 0; i < nodes[v].N_O_SZ; ++i) {
              int u = nodes[v].N_O[i];
              if (--indeg[u] == 0) {
                  q.push(u);
              }
          }
      }

    // 2. 构建树覆盖 - 优化版本
    vector<int> parent_in_tree(n, -1);
    vector<int> topo_rank(n); // 预计算拓扑序位置
    for (int i = 0; i < n; ++i) {
        topo_rank[topo_order[i]] = i;
    }
    
    for (int i = 0; i < n; ++i) {
        int v = topo_order[i];
        int max_pre = -1, max_rank = -1;
        for (int j = 0; j < nodes[v].N_I_SZ; ++j) {
            int u = nodes[v].N_I[j];
            if (u == v) continue;
            int rank = topo_rank[u]; 
            if (rank > max_rank) {
                max_rank = rank;
                max_pre = u;
            }
        }
        parent_in_tree[v] = max_pre;
    }

    
    // 3. DFS遍历 - 为当前组的节点分配时间戳
    vector<int> roots;
    roots.reserve(n/10); // 预估根节点数量
    for (int i = 0; i < n; ++i) {
        if (parent_in_tree[i] == -1) {
            roots.push_back(i);
        }
    }

    vector<int> mid(n, -1);
    vector<int> next_child(n, 0);
    int id_counter = 0;
    
    // 使用栈进行DFS遍历
    stack<int> dfs_stack;
    for (int root : roots) {
        dfs_stack.push(root);
    }
    
    while (!dfs_stack.empty()) {
        int v = dfs_stack.top();
        
        if (mid[v] != -1 && next_child[v] >= nodes[v].N_O_SZ) {
            dfs_stack.pop();
            
            
            // 只有group节点才分配新的时间戳
            if (nodes[v].group_id == group) {
                nodes[v].tout = id_counter++;
                nodes[v].t = nodes[v].tout;
                
                nodes[v].interval_array[group].reserve(1);
                nodes[v].interval_array[group].emplace_back(mid[v], nodes[v].tout);
            }
        } else {
            if (mid[v] == -1) {
                // 只有group节点才分配进入时间戳
                if (nodes[v].group_id == group) {
                    mid[v] = id_counter;
                } else {
                    mid[v] = 0; // 标记为已访问，但不分配真正的时间戳
                }
            }
            
            if (next_child[v] < nodes[v].N_O_SZ) {
                int child = nodes[v].N_O[next_child[v]++];
                if (parent_in_tree[child] == v) {
                    dfs_stack.push(child);
                }
            }
        }
    }
    
    // 使用反向拓扑序为非group节点分配时间戳（其出邻居的最小时间戳）
    for (int i = topo_order.size() - 1; i >= 0; --i) {
        int v = topo_order[i];
        
        if (nodes[v].group_id != group) {
            int min_timestamp = INT_MAX;
            
            // 遍历所有出邻居，找最小的时间戳
            for (int j = 0; j < nodes[v].N_O_SZ; ++j) {
                int neighbor = nodes[v].N_O[j];
                if (nodes[neighbor].tout < min_timestamp) {
                    min_timestamp = nodes[neighbor].tout;
                }
            }
            
            // 设置时间戳
            nodes[v].tout = min_timestamp;  // 如果没有出邻居，min_timestamp保持为INT_MAX
            
            // 非group节点不创建区间
            nodes[v].interval_array[group].clear();
        }
    }
    // 4. 区间合并 - 为每个节点构建对所有其他组的可达区间
    reverse(topo_order.begin(), topo_order.end());
    
    for (int v : topo_order) {
        if (nodes[v].N_O_SZ == 0) continue;  // 跳过叶节点
        
        // 处理每个出边
        for (int i = 0; i < nodes[v].N_O_SZ; ++i) {
            int u = nodes[v].N_O[i];
            
            if (nodes[v].group_id == group && nodes[u].group_id == group) {
                // group node到group node的前向边判断
                if (!(parent_in_tree[u] == v) && 
                    (mid[v] <= nodes[u].tout) && 
                    (nodes[u].tout < nodes[v].tout)) {
                    continue;
                }
            }
            
            // 合并区间
            if (!nodes[u].interval_array[group].empty()) {
                merge_intervals(nodes[v].interval_array[group], nodes[u].interval_array[group]);
            }
        }
    }
    
    // 清理叶节点的区间
    for (int i = 0; i < n; ++i) {
        if (nodes[i].N_O_SZ == 0) {
            
            nodes[i].interval_array[group].clear();
        
        }
    }
    
    
    // 6. 统计信息
    int total_intervals = 0;
    int max_intervals = 0;
    int nodes_with_intervals = 0;  // 新增：统计有区间的节点数

    for (int u = 0; u < n; u++) {
        int node_intervals = nodes[u].interval_array[group].size();
        
        if (node_intervals > 0) {
            nodes_with_intervals++;  // 只统计有区间的节点
        }
        
        total_intervals += node_intervals;
        max_intervals = max(max_intervals, node_intervals);
    }

    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                + (end_at.tv_usec - start_at.tv_usec) * 1e-6;

    printf("  组 %d 索引构建完成\n", group);
    printf("  - 节点数: %d\n", (int)node_groups[group].size());
    printf("  - 总区间数: %d\n", total_intervals);
    printf("  - 有区间的节点数: %d\n", nodes_with_intervals);
    printf("  - 有区间的节点平均区间数: %.2f\n", 
        nodes_with_intervals > 0 ? (double)total_intervals / nodes_with_intervals : 0.0);
    printf("  - 最大区间数: %d\n", max_intervals);
    printf("  - 构建时间: %.3f ms\n", elapsed * 1000);
}

// 构建所有组的子图（主函数）
void build_all_group_subgraphs() {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);
    int n = nodes.size();
    
    
    // 为每个组计算可达性并构建子图
    printf("\n开始逐组处理...\n");
    for (int group = 0; group < subgraph_number; group++) {
        printf("\n--- 处理组 %d/%d ---\n", group + 1, subgraph_number);
        
        
        build_ferrari_index_for_subgraph(group);
    }
    
    // 清理叶节点的区间
    for (int i = 0; i < n; ++i) {
        if (nodes[i].N_O_SZ == 0) {
            // 叶节点清空所有区间
            for (int g = 0; g < subgraph_number; ++g) {
                nodes[i].interval_array[g].clear();
            }
        }
    }
    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;
    
    // 总体统计
    printf("\n=== 组子图构建总结 ===\n");
    printf("总构建时间: %.3f ms\n", elapsed * 1000);
    
    
    
    printf("=== 组子图构建完成 ===\n\n");
}

void calculate_index_size() {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);
    
    printf("\n=== 计算索引大小和统计信息 ===\n");
    
    // 1. 计算内存占用
    size_t total_memory = 0;
    size_t interval_memory = 0;
    size_t timestamp_memory = 0;
    
    // 统计信息
    size_t total_intervals = 0;
    size_t max_intervals_per_node = 0;
    
    for (size_t i = 0; i < nodes.size(); i++) {
        size_t node_intervals = 0;
        for (int g = 0; g < subgraph_number; ++g) {
            interval_memory += nodes[i].interval_array[g].size() * sizeof(pair<int, int>);
            node_intervals += nodes[i].interval_array[g].size();
        }

        timestamp_memory += sizeof(int);  // 每个节点一个 timestamp（仍是 int）

        total_intervals += node_intervals;
        max_intervals_per_node = max(max_intervals_per_node, node_intervals);
    }
    
    total_memory = interval_memory + timestamp_memory;
    
    // 2. 输出结果
    printf("\n内存占用:\n");
    printf("  - interval_memory: %.3f MB\n", interval_memory / (1024.0 * 1024.0));
    printf("  - timestamp_memory: %.3f MB\n", timestamp_memory / (1024.0 * 1024.0));
    printf("  - total_memory: %.3f MB\n", total_memory / (1024.0 * 1024.0));
    
    printf("\n统计信息:\n");
    printf("  - total_intervals: %zu\n", total_intervals);
    printf("  - max_intervals_per_node: %zu\n", max_intervals_per_node);
    printf("  - 平均每节点区间数: %.2f\n", 
           nodes.empty() ? 0.0 : (double)total_intervals / nodes.size());
    
    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;
    printf("\n计算时间: %.3f ms\n", elapsed * 1000);
}



void verify_ferrari_index() {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);
    
    printf("\n=== 开始验证Ferrari索引 ===\n");
    
    int total_nodes = nodes.size();
    int error_count = 0;
    int total_checks = 0;
    
    
    for (int node_id = 0; node_id < total_nodes; node_id++) {
        
        vector<bool> reachable(total_nodes, false);
        stack<int> dfs_stack;
        dfs_stack.push(node_id);
        reachable[node_id] = true;
        
        while (!dfs_stack.empty()) {
            int curr = dfs_stack.top();
            dfs_stack.pop();
            
            
            for (int i = 0; i < nodes[curr].N_O_SZ; i++) {
                int next = nodes[curr].N_O[i];
                if (!reachable[next]) {
                    reachable[next] = true;
                    dfs_stack.push(next);
                }
            }
        }
        
        
        vector<vector<int>> reachable_by_group(subgraph_number);
        for (int i = 0; i < total_nodes; i++) {
            if (reachable[i] && i != node_id) {  
                int target_group = nodes[i].group_id;
                if (target_group >= 0 && target_group < subgraph_number) {
                    reachable_by_group[target_group].push_back(i);
                }
            }
        }
        
        
        for (int target_group = 0; target_group < subgraph_number; target_group++) {
            if (target_group == nodes[node_id].group_id) {
                
                continue;
            }
            
            const vector<pair<int, int>>& intervals = nodes[node_id].interval_array[target_group];
            
            bool has_error = false;
            for (int reachable_node : reachable_by_group[target_group]) {
                int tout = nodes[reachable_node].t;
                
                bool found_in_interval = false;
                for (const auto& interval : intervals) {
                    if (tout >= interval.first && tout <= interval.second) {
                        found_in_interval = true;
                        break;
                    }
                }
                
                if (!found_in_interval) {
                    if (!has_error) {
                        printf("\n节点 %d (组 %d) 的错误:\n", 
                               node_id, nodes[node_id].group_id);
                        has_error = true;
                    }
                    printf("  - 可达节点 %d (组 %d,  tout=%d) 不在任何区间内\n",
                           reachable_node, target_group, tout);
                    error_count++;
                }
                total_checks++;
            }
            
            for (const auto& interval : intervals) {
                for (int check_node : node_groups[target_group]) {
                    int tin = nodes[check_node].tin;
                    int tout = nodes[check_node].tout;
                    
                    if (tin >= interval.first && tout <= interval.second) {
                        if (!reachable[check_node]) {
                            if (!has_error) {
                                printf("\n节点 %d (组 %d) 的错误:\n", 
                                       node_id, nodes[node_id].group_id);
                                has_error = true;
                            }
                            printf("  - 不可达节点 %d (组 %d, tin=%d, tout=%d) 在区间 [%d, %d] 内\n",
                                   check_node, target_group, tin, tout, 
                                   interval.first, interval.second);
                            error_count++;
                        }
                    }
                }
            }
        }
        
        if ((node_id + 1) % 1000 == 0) {
            printf("已验证 %d/%d 个节点...\n", node_id + 1, total_nodes);
        }
    }
    
    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;
    
    // 输出验证结果
    printf("\n=== 验证完成 ===\n");
    printf("总节点数: %d\n", total_nodes);
    printf("总检查数: %d\n", total_checks);
    printf("错误数: %d\n", error_count);
    if (error_count == 0) {
        printf("✓ Ferrari索引验证通过！所有区间都正确。\n");
    } else {
        printf("✗ Ferrari索引验证失败！发现 %d 个错误。\n", error_count);
    }
    printf("验证时间: %.3f ms\n", elapsed * 1000);
    
    // 详细统计
    printf("\n详细统计:\n");
    for (int group = 0; group < subgraph_number; group++) {
        int group_total_intervals = 0;
        int group_max_intervals = 0;
        double group_avg_intervals = 0;
        
        for (int node_id : node_groups[group]) {
            int node_intervals = 0;
            for (int g = 0; g < subgraph_number; g++) {
                node_intervals += nodes[node_id].interval_array[g].size();
            }
            group_total_intervals += node_intervals;
            group_max_intervals = max(group_max_intervals, node_intervals);
        }
        
        if (!node_groups[group].empty()) {
            group_avg_intervals = (double)group_total_intervals / node_groups[group].size();
        }
        
        printf("  组 %d: 节点数=%d, 总区间数=%d, 平均区间数=%.2f, 最大区间数=%d\n",
               group, (int)node_groups[group].size(), 
               group_total_intervals, group_avg_intervals, group_max_intervals);
    }
}

void verify_node_reachability(int node_id) {
    printf("\n=== 验证节点 %d (组 %d) 的可达性 ===\n", 
           node_id, nodes[node_id].group_id);
    
    vector<bool> reachable(nodes.size(), false);
    stack<int> dfs_stack;
    dfs_stack.push(node_id);
    reachable[node_id] = true;
    
    while (!dfs_stack.empty()) {
        int curr = dfs_stack.top();
        dfs_stack.pop();
        
        for (int i = 0; i < nodes[curr].N_O_SZ; i++) {
            int next = nodes[curr].N_O[i];
            if (!reachable[next]) {
                reachable[next] = true;
                dfs_stack.push(next);
            }
        }
    }
    
    for (int target_group = 0; target_group < subgraph_number; target_group++) {
        printf("\n目标组 %d:\n", target_group);
        
        printf("  区间: ");
        const vector<pair<int, int>>& intervals = nodes[node_id].interval_array[target_group];
        if (intervals.empty()) {
            printf("无");
        } else {
            for (const auto& interval : intervals) {
                printf("[%d, %d] ", interval.first, interval.second);
            }
        }
        printf("\n");
        
        printf("  可达节点: ");
        int count = 0;
        for (int i : node_groups[target_group]) {
            if (reachable[i] && i != node_id) {
                printf("%d(tin=%d,tout=%d) ", i, nodes[i].tin, nodes[i].tout);
                count++;
                if (count >= 10) {
                    printf("...");
                    break;
                }
            }
        }
        if (count == 0) {
            printf("无");
        }
        printf("\n");
    }
}









vector<pair<node, node>> queries;

void read_queries(const char *filename) {
  queries.clear();
  timeval start_at, end_at;
  gettimeofday(&start_at, 0);
  FILE *file = fopen(filename, "r");
  int u, v;
  while (fscanf(file, "%d%d", &u, &v) == 2) {
    queries.push_back(make_pair(nodes[u], nodes[v]));
  }
  fclose(file);
  gettimeofday(&end_at, 0);
}





inline bool fast_in_intervals(const vector<pair<int,int>>& intervals, int key) {
    int n = (int)intervals.size();
    const auto* iv = intervals.data();
    if (n < BS_THRESHOLD) {
        // 线性扫描
        for (int i = 0; i < n; i++) {
            if (key >= iv[i].first && key <= iv[i].second) return true;
        }
        return false;
    }
    // 二分查找
    int l = 0, r = n - 1;
    while (l <= r) {
        int m = (l + r) >> 1;
        int lo = iv[m].first, hi = iv[m].second;
        if (key < lo)      r = m - 1;
        else if (key > hi) l = m + 1;
        else               return true;
    }
    return false;
}


inline bool reach(const node &u, const node &v) {
    if (u.id == v.id) return true;      
    int g = v.group_id;
    return fast_in_intervals(u.interval_array[g], v.t);
}


void run_queries() {
  timeval start_at, end_at;
  gettimeofday(&start_at, 0);
  int count = 0;
  for (vector<pair<node, node>>::iterator it = queries.begin(); it != queries.end(); ++it) {
    vis_cur++;
    int result = reach(it->first, it->second);
    
    if (result) {
      count++;
    }
  } 
  
  gettimeofday(&end_at, 0);
  printf("query time: %.3fms\n",
         (end_at.tv_sec - start_at.tv_sec) * 1000 +
             double(end_at.tv_usec - start_at.tv_usec) / 1000);
  printf("reachable: %d\n", count);
}

void free_all_memory() {
    // 释放每个 node 中的 N_O 和 N_I 数组
    for (auto& node : nodes) {
        delete[] node.N_O;
        delete[] node.N_I;
        node.N_O = nullptr;
        node.N_I = nullptr;

        // interval_array 是 vector 的 vector，会自动析构释放内部 pair
        node.interval_array.clear();
        node.interval_array.shrink_to_fit();
    }

    // 释放全局容器
    nodes.clear();
    nodes.shrink_to_fit();

    node_groups.clear();
    node_groups.shrink_to_fit();

    group_sizes.clear();
    group_sizes.shrink_to_fit();

    queries.clear();
    queries.shrink_to_fit();

    
}

}

int main(int argc, char *argv[]) {
  using namespace bs;
  namespace fs = std::filesystem;

  size_t last_slash = string(argv[1]).find_last_of('/');
  string parent_dir = (last_slash != string::npos) ? string(argv[1]).substr(0, last_slash) : ".";
  string dir = parent_dir + "/query";
  string prefix = "query";
  subgraph_number = stoi(argv[2]);
  string outname = parent_dir + "/result_ptil_" + to_string(subgraph_number) + "_output.log";
  if (!freopen(outname.c_str(), "w", stdout)) {
    std::fprintf(stderr, "Error: cannot redirect stdout to %s\n", outname.c_str());
    return 1;
  }
  
  read_graph(argv[1]);

	
	

	calculate_subgraph_sizes(); 
	build_all_group_subgraphs();  // 构建基于组的子图

    calculate_index_size();
    // verify_ferrari_index();
	

    for (const auto& entry : fs::directory_iterator(dir)) {
    if (!entry.is_regular_file()) continue;

    string filename = entry.path().filename().string();
    bool has_prefix = (filename.rfind(prefix, 0) == 0);

    bool ends_with_info = false;
    if (filename.size() >= 8) {
      ends_with_info = (filename.substr(filename.size() - 8) == "info.txt");
    }

    if (has_prefix && !ends_with_info) {
      string full_path = entry.path().string();
      string filename = entry.path().filename().string();
      cout << "处理文件: " << filename << endl;

      read_queries(full_path.c_str());
      run_queries();
    }
  }

  free_all_memory();
  return 0;
}