for ptil.cpp:

Compile:


g++ -O2 -std=c++17 -o ptil ptil.cpp

Run:

./ptil <graph.txt> <k.>

k. is the number of groups

query file:<graph.txt>/query/query_*.txt

for dy2.cpp:

Compile:

g++ -O2 -std=c++17 -o dy2 dy2.cpp

Run:

./dy2 <graph.txt> <per_group.> <rounds.> <K.> 

per_group. is	Initial number of nodes in each group (Vs and Vt)

rounds. is	Number of dynamic update rounds

K.	is Number of nodes randomly deleted and then inserted back per round



<graph.txt>:
 ```
# optional comment lines start with '#'
N M
u0 v01 v02 v03 ...
u1 v11 v12 ...
...
u(N-1) ...
```
- First non-comment line: two integers `N` (vertices) and `M` (edges).
- Then **N lines** follow. Line `i` starts with `i` and lists all **out-neighbors** of `i`.
- **0-based IDs**: all `u, v ∈ [0, N-1]`.
- If a vertex has no outgoing edges, write the single `u` only (no neighbors).

- **Example**
```
# toy
5 6
0 1 2
1 3
2 3 4
3 4
4
```
This encodes edges: `0→1, 0→2, 1→3, 2→3, 2→4, 3→4`.

Please note vertices are numbered from 0 to V - 1, and the graph must be a DAG.
