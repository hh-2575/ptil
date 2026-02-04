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

vector<int> group0_init_nodes;
vector<int> group1_init_nodes;

void read_graph(const char *filename) {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);

    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        return;
    }

    char buffer[1024];
    if (!fgets(buffer, sizeof(buffer), file)) {
        fprintf(stderr, "Cannot read the comment line\n");
        fclose(file);
        return;
    }

    int n, m;
    if (fscanf(file, "%d %d", &n, &m) != 2) {
        fprintf(stderr, "Failed to read n, m\n");
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

void assign_random_into_two_groups(int per_group, uint32_t seed = 0) {
    const int G = 2;
    const int total_pick = 2 * per_group;

    int n = (int)nodes.size();
    if (n < total_pick) {
        fprintf(stderr, "Error: nodes.size()=%d < %d, not enough nodes.\n", n, total_pick);
        exit(1);
    }

    node_groups.assign(G, {});
    node_groups[0].reserve(per_group);
    node_groups[1].reserve(per_group);

    group_sizes.assign(G, per_group);

    for (int u = 0; u < n; ++u) nodes[u].group_id = -1;

    vector<int> ids(n);
    iota(ids.begin(), ids.end(), 0);

    mt19937 rng;
    if (seed != 0) rng.seed(seed);
    else rng.seed(random_device{}());

    shuffle(ids.begin(), ids.end(), rng);

    for (int i = 0; i < total_pick; ++i) {
        int u = ids[i];
        int g = (i < per_group) ? 0 : 1;
        nodes[u].group_id = g;
        node_groups[g].push_back(u);
    }

    group0_init_nodes = node_groups[0];
    group1_init_nodes = node_groups[1];
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

void build_ferrari_index_Vs_to_Vt(int source_group, int target_group) {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);

    printf("Building Ferrari index for Vs(group=%d) -> Vt(group=%d)\n", source_group, target_group);

    int n = (int)nodes.size();

    for (int u = 0; u < n; ++u) {
        nodes[u].tin = -1;
        nodes[u].tout = -1;
        nodes[u].t = -1;
        if ((int)nodes[u].interval_array.size() <= target_group) {
            nodes[u].interval_array.resize(subgraph_number);
        }
        nodes[u].interval_array[target_group].clear();
    }

    vector<int> indeg(n, 0), topo_order;
    topo_order.reserve(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nodes[i].N_O_SZ; ++j) {
            indeg[nodes[i].N_O[j]]++;
        }
    }

    queue<int> q;
    for (int i = 0; i < n; ++i) if (indeg[i] == 0) q.push(i);

    while (!q.empty()) {
        int v = q.front(); q.pop();
        topo_order.push_back(v);
        for (int i = 0; i < nodes[v].N_O_SZ; ++i) {
            int u = nodes[v].N_O[i];
            if (--indeg[u] == 0) q.push(u);
        }
    }

    if ((int)topo_order.size() != n) {
        fprintf(stderr, "Error: topo_order.size()=%zu != n=%d, graph is not a DAG.\n",
                topo_order.size(), n);
        std::exit(1);
    }

    vector<int> parent_in_tree(n, -1);
    vector<int> topo_rank(n);
    for (int i = 0; i < n; ++i) topo_rank[topo_order[i]] = i;

    for (int i = 0; i < n; ++i) {
        int v = topo_order[i];
        int max_pre = -1, max_rank = -1;
        for (int j = 0; j < nodes[v].N_I_SZ; ++j) {
            int u = nodes[v].N_I[j];
            if (u == v) continue;
            int r = topo_rank[u];
            if (r > max_rank) { max_rank = r; max_pre = u; }
        }
        parent_in_tree[v] = max_pre;
    }

    vector<int> roots;
    roots.reserve(n / 10);
    for (int i = 0; i < n; ++i) if (parent_in_tree[i] == -1) roots.push_back(i);

    vector<int> mid(n, -1);
    vector<int> next_child(n, 0);
    int id_counter = 0;

    stack<int> dfs_stack;
    for (int r : roots) dfs_stack.push(r);

    while (!dfs_stack.empty()) {
        int v = dfs_stack.top();

        if (mid[v] != -1 && next_child[v] >= nodes[v].N_O_SZ) {
            dfs_stack.pop();

            if (nodes[v].group_id == target_group) {
                nodes[v].tout = id_counter++;
                nodes[v].t = nodes[v].tout;
                nodes[v].interval_array[target_group].reserve(1);
                nodes[v].interval_array[target_group].emplace_back(mid[v], nodes[v].tout);
            }
        } else {
            if (mid[v] == -1) {
                if (nodes[v].group_id == target_group) {
                    mid[v] = id_counter;
                } else {
                    mid[v] = 0;
                }
            }

            if (next_child[v] < nodes[v].N_O_SZ) {
                int child = nodes[v].N_O[next_child[v]++];
                if (parent_in_tree[child] == v) dfs_stack.push(child);
            }
        }
    }

    for (int i = (int)topo_order.size() - 1; i >= 0; --i) {
        int v = topo_order[i];
        if (nodes[v].group_id != target_group) {
            int min_timestamp = INT_MAX;
            for (int j = 0; j < nodes[v].N_O_SZ; ++j) {
                int nei = nodes[v].N_O[j];
                if (nodes[nei].tout < min_timestamp) min_timestamp = nodes[nei].tout;
            }
            nodes[v].tout = min_timestamp;
            nodes[v].interval_array[target_group].clear();
        }
    }

    reverse(topo_order.begin(), topo_order.end());

    for (int v : topo_order) {
        if (nodes[v].N_O_SZ == 0) continue;

        for (int i = 0; i < nodes[v].N_O_SZ; ++i) {
            int u = nodes[v].N_O[i];

            if (nodes[v].group_id == target_group && nodes[u].group_id == target_group) {
                if (!(parent_in_tree[u] == v) &&
                    (mid[v] <= nodes[u].tout) &&
                    (nodes[u].tout < nodes[v].tout)) {
                    continue;
                }
            }

            if (!nodes[u].interval_array[target_group].empty()) {
                merge_intervals(nodes[v].interval_array[target_group],
                                nodes[u].interval_array[target_group]);
            }
        }
    }

    for (int u = 0; u < n; ++u) {
        if (nodes[u].group_id != source_group) {
            nodes[u].interval_array[target_group].clear();
        }
    }

    int total_intervals = 0, max_intervals = 0, nodes_with_intervals = 0;
    for (int u : node_groups[source_group]) {
        int cnt = (int)nodes[u].interval_array[target_group].size();
        if (cnt > 0) nodes_with_intervals++;
        total_intervals += cnt;
        max_intervals = std::max(max_intervals, cnt);
    }

    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                  + (end_at.tv_usec - start_at.tv_usec) * 1e-6;

    printf("  Vs->Vt index construction finished\n");
    printf("  - #Vs nodes: %d\n", (int)node_groups[source_group].size());
    printf("  - Total intervals on Vs: %d\n", total_intervals);
    printf("  - #Vs nodes with intervals: %d\n", nodes_with_intervals);
    printf("  - Avg intervals (only counted nodes with intervals): %.2f\n",
           nodes_with_intervals ? (double)total_intervals / nodes_with_intervals : 0.0);
    printf("  - Max intervals: %d\n", max_intervals);
    printf("  - Construction time: %.3f ms\n", elapsed * 1000);
}

inline void remove_node_from_source_group(int source_group, int target_group, int node_id) {
    auto &vec = node_groups[source_group];
    for (int i = 0; i < (int)vec.size(); ++i) {
        if (vec[i] == node_id) {
            vec[i] = vec.back();
            vec.pop_back();
            break;
        }
    }

    group_sizes[source_group] = (int)vec.size();

    nodes[node_id].group_id = -1;

    nodes[node_id].interval_array[target_group].clear();
}

inline void add_node_to_source_group_build_label_pruned(int source_group, int target_group, int node_id) {
    nodes[node_id].group_id = source_group;
    node_groups[source_group].push_back(node_id);
    group_sizes[source_group] = (int)node_groups[source_group].size();

    auto &label = nodes[node_id].interval_array[target_group];
    label.clear();

    static vector<int> vis;
    static int vis_token = 0;
    if ((int)vis.size() != (int)nodes.size()) vis.assign(nodes.size(), 0);
    ++vis_token;

    queue<int> q;
    q.push(node_id);
    vis[node_id] = vis_token;

    vector<pair<int,int>> raw;
    raw.reserve(256);

    while (!q.empty()) {
        int v = q.front();
        q.pop();

        if (v != node_id && nodes[v].group_id == source_group) {
            const auto &vlabel = nodes[v].interval_array[target_group];
            if (!vlabel.empty()) {
                raw.insert(raw.end(), vlabel.begin(), vlabel.end());
            }
            continue;
        }

        if (nodes[v].group_id == target_group) {
            int ts = nodes[v].tout;
            if (ts != INT_MAX && ts >= 0) raw.emplace_back(ts, ts);
        }

        for (int i = 0; i < nodes[v].N_O_SZ; ++i) {
            int u = nodes[v].N_O[i];
            if (vis[u] == vis_token) continue;
            vis[u] = vis_token;
            q.push(u);
        }
    }

    if (raw.empty()) return;

    sort(raw.begin(), raw.end());

    label.reserve(raw.size());
    int L = raw[0].first, R = raw[0].second;
    for (int i = 1; i < (int)raw.size(); ++i) {
        int a = raw[i].first, b = raw[i].second;
        if (b <= R) {
            continue;
        } else if (a <= R + 1) {
            R = b;
        } else {
            label.emplace_back(L, R);
            L = a; R = b;
        }
    }
    label.emplace_back(L, R);
}

static inline void shift_intervals_after_removing_ts(vector<pair<int,int>> &intervals, int ts) {
    if (intervals.empty()) return;

    for (auto &p : intervals) {
        int l = p.first;
        int r = p.second;

        if (r < ts) {
        } else if (l >= ts) {
            p.first  = l - 1;
            p.second = r - 1;
        } else {
            p.second = r - 1;
        }
    }

    int write = 0;
    for (int i = 1; i < (int)intervals.size(); ++i) {
        if (intervals[write].second >= intervals[i].second) {
            continue;
        } else if (intervals[write].second >= intervals[i].first - 1) {
            intervals[write].second = intervals[i].second;
        } else {
            intervals[++write] = intervals[i];
        }
    }
    intervals.resize(write + 1);
}

inline void remove_node_from_target_group_shift_ts(
    int source_group,
    int target_group,
    int node_id
    ) {
    const int removed_ts = nodes[node_id].tout;

    auto &tvec = node_groups[target_group];
    for (int i = 0; i < (int)tvec.size(); ++i) {
        if (tvec[i] == node_id) {
            tvec[i] = tvec.back();
            tvec.pop_back();
            break;
        }
    }
    group_sizes[target_group] = (int)tvec.size();

    nodes[node_id].group_id = -1;
    nodes[node_id].tout = INT_MAX;
    nodes[node_id].t = -1;
    nodes[node_id].interval_array[target_group].clear();

    for (int t : node_groups[target_group]) {
        if (nodes[t].tout > removed_ts && nodes[t].tout != INT_MAX) {
            nodes[t].tout -= 1;
            nodes[t].t = nodes[t].tout;
        }
    }

    for (int s : node_groups[source_group]) {
        auto &label = nodes[s].interval_array[target_group];
        if (!label.empty()) {
            shift_intervals_after_removing_ts(label, removed_ts);
        }
    }
}

static inline void shift_intervals_after_inserting_ts(vector<pair<int,int>> &intervals, int ts) {
    if (intervals.empty()) return;

    for (auto &p : intervals) {
        int l = p.first;
        int r = p.second;

        if (r < ts) {
        } else if (l >= ts) {
            p.first  = l + 1;
            p.second = r + 1;
        } else {
            p.second = r + 1;
        }
    }

    int write = 0;
    for (int i = 1; i < (int)intervals.size(); ++i) {
        if (intervals[write].second >= intervals[i].second) {
            continue;
        } else if (intervals[write].second >= intervals[i].first - 1) {
            intervals[write].second = intervals[i].second;
        } else {
            intervals[++write] = intervals[i];
        }
    }
    intervals.resize(write + 1);
}

static inline void insert_point_into_intervals(vector<pair<int,int>> &intervals, int x) {
    if (intervals.empty()) {
        intervals.emplace_back(x, x);
        return;
    }

    auto it = std::upper_bound(
        intervals.begin(), intervals.end(), x,
        [](int val, const pair<int,int>& p){ return val < p.first; }
    );
    int pos = (int)(it - intervals.begin());

    if (pos > 0) {
        auto &L = intervals[pos - 1];
        if (L.first <= x && x <= L.second) return;
        if (L.second + 1 == x) {
            L.second = x;

            if (pos < (int)intervals.size() && intervals[pos].first == x + 1) {
                L.second = intervals[pos].second;
                intervals.erase(intervals.begin() + pos);
            }
            return;
        }
    }

    if (pos < (int)intervals.size()) {
        auto &R = intervals[pos];
        if (x + 1 == R.first) {
            R.first = x;
            return;
        }
    }

    intervals.insert(intervals.begin() + pos, {x, x});
}

inline void add_node_to_target_group_shift_ts_and_update_sources(
    int source_group,
    int target_group,
    int node_id
    ) {
    static vector<int> vis_f;
    static int vis_f_token = 0;
    if ((int)vis_f.size() != (int)nodes.size()) vis_f.assign(nodes.size(), 0);
    ++vis_f_token;

    int new_ts = INT_MAX;

    queue<int> q;
    q.push(node_id);
    vis_f[node_id] = vis_f_token;

    while (!q.empty()) {
        int v = q.front(); q.pop();

        if (v != node_id && nodes[v].group_id == target_group) {
            int ts = nodes[v].tout;
            if (ts != INT_MAX && ts >= 0) new_ts = std::min(new_ts, ts);
        }

        for (int i = 0; i < nodes[v].N_O_SZ; ++i) {
            int u = nodes[v].N_O[i];
            if (vis_f[u] == vis_f_token) continue;
            vis_f[u] = vis_f_token;
            q.push(u);
        }
    }

    if (new_ts == INT_MAX) {
        new_ts = (int)node_groups[target_group].size();
    }

    for (int t : node_groups[target_group]) {
        if (nodes[t].tout != INT_MAX && nodes[t].tout >= new_ts) {
            nodes[t].tout += 1;
            nodes[t].t = nodes[t].tout;
        }
    }

    nodes[node_id].group_id = target_group;
    nodes[node_id].tout = new_ts;
    nodes[node_id].t = new_ts;
    nodes[node_id].interval_array[target_group].clear();

    node_groups[target_group].push_back(node_id);
    group_sizes[target_group] = (int)node_groups[target_group].size();

    for (int s : node_groups[source_group]) {
        auto &label = nodes[s].interval_array[target_group];
        if (!label.empty()) shift_intervals_after_inserting_ts(label, new_ts);
    }

    static vector<int> vis_b;
    static int vis_b_token = 0;
    if ((int)vis_b.size() != (int)nodes.size()) vis_b.assign(nodes.size(), 0);
    ++vis_b_token;

    queue<int> rq;
    rq.push(node_id);
    vis_b[node_id] = vis_b_token;

    while (!rq.empty()) {
        int v = rq.front(); rq.pop();

        if (nodes[v].group_id == source_group) {
            auto &label = nodes[v].interval_array[target_group];
            insert_point_into_intervals(label, new_ts);
        }

        for (int i = 0; i < nodes[v].N_I_SZ; ++i) {
            int u = nodes[v].N_I[i];
            if (vis_b[u] == vis_b_token) continue;
            vis_b[u] = vis_b_token;
            rq.push(u);
        }
    }
}

static inline long long now_us() {
    timeval tv; gettimeofday(&tv, nullptr);
    return (long long)tv.tv_sec * 1000000LL + (long long)tv.tv_usec;
}

static inline vector<int> sample_k_unassigned(int k, mt19937 &rng) {
    vector<int> cand;
    cand.reserve(nodes.size());
    for (int u = 0; u < (int)nodes.size(); ++u) {
        if (nodes[u].group_id == -1) cand.push_back(u);
    }
    if ((int)cand.size() < k) {
        fprintf(stderr, "[ERR] not enough unassigned nodes: have %zu, need %d\n",
                cand.size(), k);
        exit(1);
    }
    shuffle(cand.begin(), cand.end(), rng);
    cand.resize(k);
    return cand;
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
}

static inline bool interval_contains(const vector<pair<int,int>> &intervals, int x) {
    if (intervals.empty()) return false;
    auto it = std::upper_bound(
        intervals.begin(), intervals.end(), x,
        [](int val, const pair<int,int>& p){ return val < p.first; }
    );
    if (it == intervals.begin()) return false;
    --it;
    return it->first <= x && x <= it->second;
}

static inline void collect_reachable_targets_from_source(
    int s,
    int target_group,
    vector<int> &out_targets,
    vector<int> &vis,
    int &vis_token
    ) {
    out_targets.clear();
    ++vis_token;

    std::stack<int> st;
    st.push(s);
    vis[s] = vis_token;

    while (!st.empty()) {
        int v = st.top();
        st.pop();

        if (nodes[v].group_id == target_group) {
            out_targets.push_back(v);
        }

        for (int i = 0; i < nodes[v].N_O_SZ; ++i) {
            int u = nodes[v].N_O[i];
            if (vis[u] == vis_token) continue;
            vis[u] = vis_token;
            st.push(u);
        }
    }
}

bool verify_Vs_to_Vt_completeness(
    int source_group,
    int target_group,
    int max_report = 20,
    bool print_progress = true
    ) {
    timeval start_at, end_at;
    gettimeofday(&start_at, nullptr);

    int n = (int)nodes.size();

    int bad_target_ts = 0;
    for (int t : node_groups[target_group]) {
        if (nodes[t].tout == INT_MAX || nodes[t].tout < 0) bad_target_ts++;
    }
    if (bad_target_ts > 0) {
        printf("[WARN] %d targets in group %d have invalid tout (INT_MAX or <0). "
               "They will break verification.\n", bad_target_ts, target_group);
    }

    vector<int> vis(n, 0);
    int vis_token = 0;

    vector<int> reachable_targets;
    long long total_pairs_checked = 0;
    int total_missing = 0;
    int reported = 0;

    const int Vs_size = (int)node_groups[source_group].size();
    for (int idx = 0; idx < Vs_size; ++idx) {
        int s = node_groups[source_group][idx];

        if (print_progress && (idx % 50 == 0)) {
            printf("[verify] %d / %d sources\n", idx, Vs_size);
        }

        collect_reachable_targets_from_source(s, target_group, reachable_targets, vis, vis_token);

        const auto &label = nodes[s].interval_array[target_group];

        for (int t : reachable_targets) {
            int ts = nodes[t].tout;
            if (ts == INT_MAX || ts < 0) {
                total_missing++;
                if (reported < max_report) {
                    printf("[FAIL] s=%d -> t=%d reachable, but t has invalid tout=%d\n", s, t, ts);
                    reported++;
                }
                continue;
            }

            total_pairs_checked++;
            if (!interval_contains(label, ts)) {
                total_missing++;
                if (reported < max_report) {
                    printf("[FAIL] Missing coverage: s=%d can reach t=%d (tout=%d), "
                           "but tout not in s.label intervals (size=%zu)\n",
                           s, t, ts, label.size());

                    int print_k = (int)std::min<size_t>(label.size(), 10);
                    printf("       s.label first %d intervals: ", print_k);
                    for (int i = 0; i < print_k; ++i) {
                        printf("[%d,%d] ", label[i].first, label[i].second);
                    }
                    printf("\n");
                    reported++;
                }
            }
        }
    }

    gettimeofday(&end_at, nullptr);
    double elapsed = (end_at.tv_sec - start_at.tv_sec)
                   + (end_at.tv_usec - start_at.tv_usec) * 1e-6;

    printf("=== Verification summary (completeness) ===\n");
    printf("Vs size: %d, Vt size: %d\n",
           (int)node_groups[source_group].size(),
           (int)node_groups[target_group].size());
    printf("Total reachable (s,t) pairs checked: %lld\n", total_pairs_checked);
    printf("Total missing coverages (false negatives): %d\n", total_missing);
    printf("Time: %.3f ms\n", elapsed * 1000);

    return total_missing == 0;
}

}

int main(int argc, char *argv[]) {

    if (argc < 5) {
        fprintf(stderr, "Usage: %s <graph_path> <per_group> <rounds> <K> [seed]\n", argv[0]);
        return 1;
    }

    int per_group = atoi(argv[2]);
    int ROUNDS    = atoi(argv[3]);
    int K         = atoi(argv[4]);
    uint32_t seed = (argc >= 6) ? (uint32_t)strtoul(argv[5], nullptr, 10) : 0;

    if (per_group <= 0 || ROUNDS <= 0 || K <= 0) {
        fprintf(stderr, "Error: per_group/rounds/K must be positive.\n");
        return 1;
    }

    using namespace bs;
    namespace fs = std::filesystem;

    size_t last_slash = string(argv[1]).find_last_of('/');
    string parent_dir = (last_slash != string::npos) ? string(argv[1]).substr(0, last_slash) : ".";
    subgraph_number = 2;
    string outname = parent_dir + "/result_dy2_"
                 + to_string(per_group) + "_"
                 + to_string(ROUNDS) + "_"
                 + to_string(K)
                 + "_output.log";
    if (!freopen(outname.c_str(), "w", stdout)) {
        std::fprintf(stderr, "Error: cannot redirect stdout to %s\n", outname.c_str());
        return 1;
    }

    read_graph(argv[1]);

    assign_random_into_two_groups(per_group, seed);

    build_ferrari_index_Vs_to_Vt(0, 1);

    mt19937 rng((uint32_t)random_device{}());

    double g0_del_total_ms = 0.0, g0_add_total_ms = 0.0;

    for (int round = 1; round <= ROUNDS; ++round) {
        int cur0 = (int)node_groups[0].size();
        int del_k = min(K, cur0);

        vector<int> to_remove = node_groups[0];
        shuffle(to_remove.begin(), to_remove.end(), rng);
        to_remove.resize(del_k);

        long long t1 = now_us();
        for (int s : to_remove) remove_node_from_source_group(0, 1, s);
        long long t2 = now_us();
        double del_ms = (t2 - t1) / 1000.0;
        g0_del_total_ms += del_ms;

        vector<int> to_add = sample_k_unassigned(del_k, rng);

        long long t3 = now_us();
        for (int s : to_add) add_node_to_source_group_build_label_pruned(0, 1, s);
        long long t4 = now_us();
        double add_ms = (t4 - t3) / 1000.0;
        g0_add_total_ms += add_ms;

        printf("[Group0] round %d: delete=%d %.3f ms, add=%d %.3f ms\n",
            round, del_k, del_ms, del_k, add_ms);
    }

    printf("=== Group0 summary (%d rounds) ===\n", ROUNDS);
    printf("Total delete time: %.3f ms\n", g0_del_total_ms);
    printf("Total add time   : %.3f ms\n", g0_add_total_ms);

    double g1_del_total_ms = 0.0, g1_add_total_ms = 0.0;

    for (int round = 1; round <= ROUNDS; ++round) {
        int cur1 = (int)node_groups[1].size();
        int del_k = min(K, cur1);

        vector<int> to_remove = node_groups[1];
        shuffle(to_remove.begin(), to_remove.end(), rng);
        to_remove.resize(del_k);

        long long t1 = now_us();
        for (int t : to_remove) remove_node_from_target_group_shift_ts(0, 1, t);
        long long t2 = now_us();
        double del_ms = (t2 - t1) / 1000.0;
        g1_del_total_ms += del_ms;

        vector<int> to_add = sample_k_unassigned(del_k, rng);

        long long t3 = now_us();
        for (int t : to_add) add_node_to_target_group_shift_ts_and_update_sources(0, 1, t);
        long long t4 = now_us();
        double add_ms = (t4 - t3) / 1000.0;
        g1_add_total_ms += add_ms;

        printf("[Group1] round %d: delete=%d %.3f ms, add=%d %.3f ms\n",
            round, del_k, del_ms, del_k, add_ms);
    }

    printf("=== Group1 summary (%d rounds) ===\n", ROUNDS);
    printf("Total delete time: %.3f ms\n", g1_del_total_ms);
    printf("Total add time   : %.3f ms\n", g1_add_total_ms);

    free_all_memory();
    return 0;
}
