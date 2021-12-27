// Max Flow by Dinic
//#pragma GCC optimize(2)
#include<iostream>
#include<vector>
#include<list>
#include<queue>

using namespace std;

template<class T>
struct edge{
    // variable id records both whether it's a reverse edge and its index in the adjaency list
    // if   id >= 0 , it is a forward edge and graph[i][id] is this edge
    // if   id <  0 , it is a reverse edge and graph[next][-id-1] is this edge  
    int id;  

    // the start of the edge is unneccessary to record, only the end is needed
    int next;
    T val;
};


// Construct a residual graph given the network graph (graph, capacity, flow)
template<class T>
inline vector<list<struct edge<T>>> ResidualGraph(const vector<vector<int>>& graph,
                                        const vector<vector<T>>& capacity,
                                        const vector<vector<T>>& flow){
    int n = graph.size();
    vector<list<struct edge<T>>> residualgraph(n);
    for (int i=0;i<n;++i){
        int m = graph[i].size();
        for (int j=0;j<m;++j){
            if (graph[i][j] == i) continue;
            if (flow[i][j] > 0){
                // reverse edge
                residualgraph[graph[i][j]].push_back((struct edge<T>){-j-1,i,flow[i][j]});
            }
            if (flow[i][j] != capacity[i][j]){
                residualgraph[i].push_back(
                        (struct edge<T>){j,graph[i][j],capacity[i][j] - flow[i][j]});
            }
        }
    }
    return residualgraph;
}


// Modify the residual graph to a level graph by erasing some of the edges
template<class T>
void LevelGraph(vector<list<struct edge<T>>>& resgraph, int source,int sink){
    typedef typename list<struct edge<T>>::iterator itertype;
    int n = resgraph.size();
    deque<int> q;
    vector<int> level(n,n+2);
    q.push_back(source);
    level[source] = 1;
    while (!q.empty()){
        int i = q.front();
        q.pop_front();
        if (i == sink){
            break;
        }
        for (struct edge<T>& neighbor: resgraph[i]){
            if (level[neighbor.next] > n){
                level[neighbor.next] = level[i] + 1;
                q.push_back(neighbor.next);
            }
        }
    }

    for (int i=0;i<n;++i){
        for (itertype iter=resgraph[i].begin(); iter!=resgraph[i].end();){
            if (level[iter->next] != level[i] + 1){
                iter = resgraph[i].erase(iter);
            }else{
                ++iter;
            }
        }
    }

    return ;
}


// Search a feasible path on the levelgraph
// It is guaranteed that there is no cycle on the (directed) levelgraph
template<class T>
bool DFS(vector<list<struct edge<T>>>& levelgraph,
                    vector<typename list<struct edge<T>>::iterator>& path,
                    int i,int sink){
    typedef typename list<struct edge<T>>::iterator itertype;
    for (itertype iter=levelgraph[i].begin(); iter!=levelgraph[i].end();){
        if (iter->val == T(0)){  // zero capacity left, erase it!
            iter = levelgraph[i].erase(iter);
        }else if (iter->next == sink){  // reach the sink (target)
            path.push_back(iter);
            return true;          
        }else if (levelgraph[iter->next].empty()){  // dead end, erase it!
            iter = levelgraph[i].erase(iter);
        }else{ // go dfs
            path.push_back(iter);
            if (DFS(levelgraph,path,iter->next,sink)) return true;
            path.pop_back();
            ++iter;
        }
    }
    return false;
}


// Search a feasible path in the levelgraph from source to sink
// and overwrite the path in the vector ``path``
template<class T>
void DFSHandler(vector<list<struct edge<T>>>& levelgraph,
                vector<typename list<struct edge<T>>::iterator>& path,
                int source,int sink){
    path.clear();
    DFS(levelgraph,path,source,sink);
}


// @brief Perform Dinic algorithm on the network, self-cycles and duplicated edges are permitted
// @param graph: Adjacency List, where graph[i] records the indices that can reach from i
// @param capacity: Same shape as graph, capacity[i][j] records the capacity of edge[i][j]
// @param source: Source index
// @param sink: Sink index
// @return the value of max flow
template<class T>
T Dinic(vector<vector<int>>& graph,
            vector<vector<T>>& capacity, int source, int sink){
    if (source == sink) return T(0);
    // Initialization
    int n = graph.size();
    vector<vector<T>> flow(n);
    vector<typename list<struct edge<T>>::iterator> path;
    for (int i=0;i<n;++i){
        flow[i] = vector<T>(graph[i].size());
    }

    bool flg = 1;
    while (flg){ // Phase
        flg = 0;
        // Create Level Graph of the Residual Graph
        vector<list<struct edge<T>>> levelgraph = ResidualGraph(graph,capacity,flow);
        LevelGraph(levelgraph,source,sink); // Convert the residual graph to level graph
        
        while (1){
            // DFS
            DFSHandler<T>(levelgraph,path,source,sink);
            if (path.empty()) break;
            flg = 1;

            T maxflow = path[0]->val;
            for (typename list<struct edge<T>>::iterator e: path){
                maxflow = min(maxflow, e->val);
            }
 
            // Augment
            int start = source;
            for (typename list<struct edge<T>>::iterator e: path){
                if (e->id >= 0){
                    // it is not a reverse edge
                    flow[start][e->id] += maxflow;

                    // update level graph
                    // (reverse edge is not in level graph in this case)
                    e->val -= maxflow;
                }else{
                    // it is a reverse edge
                    flow[e->next][-e->id - 1] -= maxflow;

                    // update level graph
                    // (forward edge is not in level graph in this case)
                    e->val += maxflow;
                }
                start = e->next;
            }
        } 
    }    
    
    // print the flow
    //for (int i=0;i<n;++i) for (int j=graph[i].size()-1;j>=0;--j) 
    //  if (flow[i][j]) cout << (i+1) << "->" << (graph[i][j]+1) << ": " << flow[i][j] <<'\n';

    T ans(0);
    for (T out: flow[source]){
        if (out) ans += out;
    }
    return ans;
}

// Luogu P3376
long long merge[205][205]; // 合并重边 (不合并亦可)
int main(){
    int n,m,s,t;
    ios::sync_with_stdio(false);
    //freopen("D:\\CppProjects\\LuoGu\\P3376_9.in","r",stdin);
    cin >> n >> m >> s >> t;
    --s, --t;
    vector<vector<int>> graph(n);
    vector<vector<long long>> capacity(n);
    int u,v;
    long long w;
    for (int i=0;i<m;++i){
        cin >> u >> v >> w;
        --u, --v;
        merge[u][v] += w;
        graph[u].push_back(v);
        capacity[u].push_back(w);
    }
    
    for (int i=0;i<n;++i){
        for (int j=0;j<n;++j){
            if (merge[i][j]){
                graph[i].push_back(j);
                capacity[i].push_back(merge[i][j]);
            }
        }
    }
    
    cout << Dinic(graph,capacity,s,t);
    return 0;
}
