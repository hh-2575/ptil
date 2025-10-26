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
  int group_id;
  vector<vector<pair<int, int>>> interval_array;
  int tin, tout;
  int t;
};

vector<node> nodes;
vector<vector<int>> node_groups;

int vis_cur, cur;
int subgraph_number;
vector<int> group_sizes;

static constexpr int BS_THRESHOLD = 8;

void read_graph(const char *filename) {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);

    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Failed to open file %s\n", filename);
        return;
    }

    char buffer[1024];
    if (!fgets(buffer, sizeof(buffer), file)) {
        fprintf(stderr, "Failed to read the comment line\n");
        fclose(file);
        return;
    }

    int n, m;
    if (fscanf(file, "%d %d", &n, &m) != 2) {
        fprintf(stderr, "Failed to read n and m\n");
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
        nodes[u].id = u;
        nodes[u].N_O_SZ = (int)N_O[u].size();
        nodes[u].N_O = new int[nodes[u].N_O_SZ];
        for (int i = 0; i < nodes[u].N_O_SZ; i++)
            nodes[u].N_O[i] = N_O[u][i];

        nodes[u].N_I_SZ = (int)N_I[u].size();
        nodes[u].N_I = new int[nodes[u].N_I_SZ];
        for (int i = 0; i < nodes[u].N_I_SZ; i++)
            nodes[u].N_I[i] = N_I[u][i];

        nodes[u].interval_array.resize(subgraph_number);
        nodes[u].vis = 0;
        nodes[u].group_id = -1;
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

    printf("\n=== Compute subgraph group sizes (random even partition of all nodes) ===\n");

    int total_nodes = nodes.size();
    printf("Total nodes: %d\n", total_nodes);
    printf("Number of groups: %d\n", subgraph_number);

    if (subgraph_number <= 0) {
        printf("Error: subgraph_number must be > 0\n");
        return;
    }

    if (total_nodes == 0) {
        printf("Error: no nodes to partition\n");
        return;
    }

    group_sizes.resize(subgraph_number);
    node_groups.resize(subgraph_number);

    int base_size = total_nodes / subgraph_number;
    int remainder = total_nodes % subgraph_number;

    for (int i = 0; i < subgraph_number; i++) {
        if (i < remainder) {
            group_sizes[i] = base_size + 1;
        } else {
            group_sizes[i] = base_size;
        }
        node_groups[i].reserve(group_sizes[i]);
    }

    vector<int> node_ids(total_nodes);
    for (int i = 0; i < total_nodes; i++) {
        node_ids[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(node_ids.begin(), node_ids.end(), gen);

    int node_index = 0;
    for (int group = 0; group < subgraph_number; group++) {
        for (int i = 0; i < group_sizes[group]; i++) {
            int node_id = node_ids[node_index++];
            nodes[node_id].group_id = group;
            node_groups[group].push_back(node_id);
        }
    }

    printf("\nGrouping result (random even split):\n");
    int total_check = 0;
    int min_size = INT_MAX;
    int max_size = 0;

    for (int i = 0; i < subgraph_number; i++) {
        printf("  Group %d: %d nodes\n", i, group_sizes[i]);
        total_check += group_sizes[i];
        min_size = min(min_size, group_sizes[i]);
        max_size = max(max_size, group_sizes[i]);
    }

    printf("\nVerify:\n");
    printf("  Total assigned nodes: %d\n", total_check);
    printf("  Expected nodes: %d\n", total_nodes);
    printf("  Min group size: %d\n", min_size);
    printf("  Max group size: %d\n", max_size);
    printf("  Size gap: %d (should be <= 1)\n", max_size - min_size);

    int unassigned_count = 0;
    for (int i = 0; i < total_nodes; i++) {
        if (nodes[i].group_id == -1) {
            unassigned_count++;
        }
    }
    if (unassigned_count > 0) {
        printf("  Warning: %d nodes unassigned!\n", unassigned_count);
    }

    double mean = (double)total_nodes / subgraph_number;
    double variance = 0;
    for (int i = 0; i < subgraph_number; i++) {
        variance += pow(group_sizes[i] - mean, 2);
    }
    variance /= subgraph_number;
    double std_dev = sqrt(variance);

    printf("  Average group size: %.2f\n", mean);
    printf("  Std dev: %.2f (smaller is more even)\n", std_dev);

    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;
    printf("\nElapsed: %.3f ms\n", elapsed * 1000);
}

void merge_intervals(vector<pair<int,int>> &to, const vector<pair<int,int>> &from) {
    if (from.empty()) return;
    if (to.empty()) {
        to = from;
        return;
    }
    to.reserve(to.size() + from.size());
    to.insert(to.end(), from.begin(), from.end());
    sort(to.begin(), to.end());
    int write_idx = 0;
    for (int read_idx = 1; read_idx < (int)to.size(); ++read_idx) {
        if (to[write_idx].second >= to[read_idx].second) {
            continue;
        } else if (to[write_idx].second >= to[read_idx].first - 1) {
            to[write_idx].second = max(to[write_idx].second, to[read_idx].second);
        } else {
            to[++write_idx] = to[read_idx];
        }
    }
    to.resize(write_idx + 1);
}

void build_ferrari_index_for_subgraph(int group) {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);

    printf("Build Ferrari index for group %d\n", group);

    int n = nodes.size();

    for (int u = 0; u < n; u++) {
        nodes[u].tin = -1;
        nodes[u].tout = -1;
    }

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

    vector<int> parent_in_tree(n, -1);
    vector<int> topo_rank(n);
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

    vector<int> roots;
    roots.reserve(n/10);
    for (int i = 0; i < n; ++i) {
        if (parent_in_tree[i] == -1) {
            roots.push_back(i);
        }
    }

    vector<int> mid(n, -1);
    vector<int> next_child(n, 0);
    int id_counter = 0;

    stack<int> dfs_stack;
    for (int root : roots) {
        dfs_stack.push(root);
    }

    while (!dfs_stack.empty()) {
        int v = dfs_stack.top();

        if (mid[v] != -1 && next_child[v] >= nodes[v].N_O_SZ) {
            dfs_stack.pop();
            if (nodes[v].group_id == group) {
                nodes[v].tout = id_counter++;
                nodes[v].t = nodes[v].tout;
                nodes[v].interval_array[group].reserve(1);
                nodes[v].interval_array[group].emplace_back(mid[v], nodes[v].tout);
            }
        } else {
            if (mid[v] == -1) {
                if (nodes[v].group_id == group) {
                    mid[v] = id_counter;
                } else {
                    mid[v] = 0;
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

    for (int i = (int)topo_order.size() - 1; i >= 0; --i) {
        int v = topo_order[i];

        if (nodes[v].group_id != group) {
            int min_timestamp = INT_MAX;
            for (int j = 0; j < nodes[v].N_O_SZ; ++j) {
                int neighbor = nodes[v].N_O[j];
                if (nodes[neighbor].tout < min_timestamp) {
                    min_timestamp = nodes[neighbor].tout;
                }
            }
            nodes[v].tout = min_timestamp;
            nodes[v].interval_array[group].clear();
        }
    }

    reverse(topo_order.begin(), topo_order.end());

    for (int v : topo_order) {
        if (nodes[v].N_O_SZ == 0) continue;
        for (int i = 0; i < nodes[v].N_O_SZ; ++i) {
            int u = nodes[v].N_O[i];

            if (nodes[v].group_id == group && nodes[u].group_id == group) {
                if (!(parent_in_tree[u] == v) &&
                    (mid[v] <= nodes[u].tout) &&
                    (nodes[u].tout < nodes[v].tout)) {
                    continue;
                }
            }

            if (!nodes[u].interval_array[group].empty()) {
                merge_intervals(nodes[v].interval_array[group], nodes[u].interval_array[group]);
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        if (nodes[i].N_O_SZ == 0) {
            nodes[i].interval_array[group].clear();
        }
    }

    int total_intervals = 0;
    int max_intervals = 0;
    int nodes_with_intervals = 0;

    for (int u = 0; u < n; u++) {
        int node_intervals = nodes[u].interval_array[group].size();
        if (node_intervals > 0) {
            nodes_with_intervals++;
        }
        total_intervals += node_intervals;
        max_intervals = max(max_intervals, node_intervals);
    }

    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                + (end_at.tv_usec - start_at.tv_usec) * 1e-6;

    printf("  Group %d index built\n", group);
    printf("  - number of nodes: %d\n", (int)node_groups[group].size());
    printf("  - total intervals: %d\n", total_intervals);
    printf("  - nodes with intervals: %d\n", nodes_with_intervals);
    printf("  - avg intervals on nodes with intervals: %.2f\n",
        nodes_with_intervals > 0 ? (double)total_intervals / nodes_with_intervals : 0.0);
    printf("  - max intervals: %d\n", max_intervals);
    printf("  - build time: %.3f ms\n", elapsed * 1000);
}

void build_all_group_subgraphs() {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);
    int n = nodes.size();

    printf("\nStart processing groups...\n");
    for (int group = 0; group < subgraph_number; group++) {
        printf("\n--- Processing group %d/%d ---\n", group + 1, subgraph_number);
        build_ferrari_index_for_subgraph(group);
    }

    for (int i = 0; i < n; ++i) {
        if (nodes[i].N_O_SZ == 0) {
            for (int g = 0; g < subgraph_number; ++g) {
                nodes[i].interval_array[g].clear();
            }
        }
    }
    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;

    printf("\n=== Summary of group subgraph construction ===\n");
    printf("Total build time: %.3f ms\n", elapsed * 1000);
    printf("=== Finished building group subgraphs ===\n\n");
}

void calculate_index_size() {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);

    printf("\n=== Compute index size and statistics ===\n");

    size_t total_memory = 0;
    size_t interval_memory = 0;
    size_t timestamp_memory = 0;

    size_t total_intervals = 0;
    size_t max_intervals_per_node = 0;

    for (size_t i = 0; i < nodes.size(); i++) {
        size_t node_intervals = 0;
        for (int g = 0; g < subgraph_number; ++g) {
            interval_memory += nodes[i].interval_array[g].size() * sizeof(pair<int, int>);
            node_intervals += nodes[i].interval_array[g].size();
        }

        timestamp_memory += sizeof(int);

        total_intervals += node_intervals;
        max_intervals_per_node = max(max_intervals_per_node, node_intervals);
    }

    total_memory = interval_memory + timestamp_memory;

    printf("\nMemory usage:\n");
    printf("  - interval_memory: %.3f MB\n", interval_memory / (1024.0 * 1024.0));
    printf("  - timestamp_memory: %.3f MB\n", timestamp_memory / (1024.0 * 1024.0));
    printf("  - total_memory: %.3f MB\n", total_memory / (1024.0 * 1024.0));

    printf("\nStatistics:\n");
    printf("  - total_intervals: %zu\n", total_intervals);
    printf("  - max_intervals_per_node: %zu\n", max_intervals_per_node);
    printf("  - avg intervals per node: %.2f\n",
           nodes.empty() ? 0.0 : (double)total_intervals / nodes.size());

    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;
    printf("\nElapsed: %.3f ms\n", elapsed * 1000);
}

void verify_ferrari_index() {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);

    printf("\n=== Start verifying Ferrari index ===\n");

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
                        printf("\nErrors for node %d (group %d):\n",
                               node_id, nodes[node_id].group_id);
                        has_error = true;
                    }
                    printf("  - reachable node %d (group %d, tout=%d) not covered by any interval\n",
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
                                printf("\nErrors for node %d (group %d):\n",
                                       node_id, nodes[node_id].group_id);
                                has_error = true;
                            }
                            printf("  - unreachable node %d (group %d, tin=%d, tout=%d) lies in interval [%d, %d]\n",
                                   check_node, target_group, tin, tout,
                                   interval.first, interval.second);
                            error_count++;
                        }
                    }
                }
            }
        }

        if ((node_id + 1) % 1000 == 0) {
            printf("Verified %d/%d nodes...\n", node_id + 1, total_nodes);
        }
    }

    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;

    printf("\n=== Verification complete ===\n");
    printf("Total nodes: %d\n", total_nodes);
    printf("Total checks: %d\n", total_checks);
    printf("Errors: %d\n", error_count);
    if (error_count == 0) {
        printf("✓ Ferrari index verified: all intervals are correct.\n");
    } else {
        printf("✗ Ferrari index verification failed: %d errors found.\n", error_count);
    }
    printf("Verification time: %.3f ms\n", elapsed * 1000);

    printf("\nDetailed stats:\n");
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

        printf("  Group %d: nodes=%d, total intervals=%d, avg intervals=%.2f, max intervals=%d\n",
               group, (int)node_groups[group].size(),
               group_total_intervals, group_avg_intervals, group_max_intervals);
    }
}

void verify_node_reachability(int node_id) {
    printf("\n=== Verify reachability for node %d (group %d) ===\n",
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
        printf("\nTarget group %d:\n", target_group);

        printf("  Intervals: ");
        const vector<pair<int, int>>& intervals = nodes[node_id].interval_array[target_group];
        if (intervals.empty()) {
            printf("None");
        } else {
            for (const auto& interval : intervals) {
                printf("[%d, %d] ", interval.first, interval.second);
            }
        }
        printf("\n");

        printf("  Reachable nodes: ");
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
            printf("None");
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
        for (int i = 0; i < n; i++) {
            if (key >= iv[i].first && key <= iv[i].second) return true;
        }
        return false;
    }
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
    for (auto& node : nodes) {
        delete[] node.N_O;
        delete[] node.N_I;
        node.N_O = nullptr;
        node.N_I = nullptr;
        node.interval_array.clear();
        node.interval_array.shrink_to_fit();
    }

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
  build_all_group_subgraphs();

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
      cout << "Processing file: " << filename << endl;

      read_queries(full_path.c_str());
      run_queries();
    }
  }

  free_all_memory();
  return 0;
}
