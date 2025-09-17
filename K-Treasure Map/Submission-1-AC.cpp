#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int n, m, k, tx, ty;
    if(!(cin >> n >> m >> k >> tx >> ty)) return 0;
    int N = n + m;
    
    vector<vector<pair<int,long long>>> adj(N+1);
    for(int i=0;i<k;i++){
        int x,y; long long d;
        cin >> x >> y >> d;
        int u = x;
        int v = n + y;
        adj[u].push_back({v,d});
        adj[v].push_back({u,d});
    }
    
    vector<int> comp(N+1,-1);
    vector<long long> val(N+1,0);
    vector<char> vis(N+1,0);
    int cid=0;
    const long long INF = (1LL<<62);
    vector<long long> minR(1), minS(1);
    vector<char> hasX(1), hasY(1);
    
    for(int u=1; u<=N; ++u){
        if(vis[u]) continue;
        ++cid;
        if((int)minR.size()<=cid){
            minR.push_back(INF);
            minS.push_back(INF);
            hasX.push_back(0);
            hasY.push_back(0);
        }
        val[u]=0; vis[u]=1; comp[u]=cid;
        deque<int> dq; dq.push_back(u);
        while(!dq.empty()){
            int x = dq.front(); dq.pop_front();
            if(x <= n){ hasX[cid]=1; minR[cid] = min(minR[cid], val[x]); }
            else { hasY[cid]=1; minS[cid] = min(minS[cid], val[x]); }
            for(auto &e: adj[x]){
                int v = e.first; long long w = e.second;
                long long need = w - val[x];
                if(!vis[v]){ vis[v]=1; comp[v]=cid; val[v]=need; dq.push_back(v); }
                else { if(val[v] != need){ cout << "impossible\n"; return 0; } }
            }
        }
    }
    
    for(int c=1;c<=cid;c++){
        if(hasX[c] && hasY[c]){
            if(minR[c] + minS[c] < 0){
                cout << "impossible\n";
                return 0;
            }
        }
    }
    
    int nodeX = tx;
    int nodeY = n + ty;
    long long base = val[nodeX] + val[nodeY];
    if(comp[nodeX] == comp[nodeY]){
        cout << base << "\n";
    } else {
        long long delta = -minR[comp[nodeX]] - minS[comp[nodeY]];
        cout << base + delta << "\n";
    }
    return 0;
}