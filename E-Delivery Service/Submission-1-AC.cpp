#include <bits/stdc++.h>
using namespace std;

struct DSU {
    int N;
    vector<int> parent;
    vector<int> sz; // number of nodes in component
    vector<long long> f; // number of cities whose both nodes in this component
    vector<unordered_map<int,int>> mp; // counts of cities split across pairs
    long long tot1 = 0, tot2 = 0;

    DSU(int n=0){ init(n); }

    void init(int n){
        N = n;
        parent.resize(N);
        sz.assign(N,1);
        f.assign(N,0);
        mp.assign(N, {});
        for(int i=0;i<N;i++){ parent[i]=i; mp[i].reserve(2); }
        tot1 = 0; tot2 = 0;
    }

    static inline long long C2(long long x){ return x*(x-1)/2; }

    int find(int x){
        while(parent[x]!=x){
            parent[x]=parent[parent[x]];
            x=parent[x];
        }
        return x;
    }

    void add_split_pair(int a, int b){
        mp[a][b]++; mp[b][a]++;
    }

    void unite(int u, int v){
        int ru=find(u), rv=find(v);
        if(ru==rv) return;
        if(mp[ru].size() + (size_t)sz[ru] < mp[rv].size() + (size_t)sz[rv]) swap(ru, rv);

        // handle pair between ru and rv
        long long between = 0;
        auto it = mp[ru].find(rv);
        if (it != mp[ru].end()){
            between = it->second;
            tot2 -= C2(between);
            mp[ru].erase(it);
            auto it2 = mp[rv].find(ru);
            if (it2 != mp[rv].end()) mp[rv].erase(it2);
        } else {
            auto it2 = mp[rv].find(ru);
            if (it2 != mp[rv].end()){
                between = it2->second;
                tot2 -= C2(between);
                mp[rv].erase(it2);
            }
        }

        // Update tot1 for components
        tot1 -= C2((long long)sz[ru]-f[ru]) + C2((long long)sz[rv]-f[rv]);

        parent[rv]=ru;
        sz[ru]+=sz[rv];
        f[ru]+=f[rv]+between;

        // Merge maps
        for(auto &kv : mp[rv]){
            int w = kv.first;
            if(w==ru) continue;
            int cnt = kv.second;

            long long old = 0;
            auto itaw = mp[ru].find(w);
            if(itaw != mp[ru].end()) old = itaw->second;

            tot2 -= C2(old) + C2(cnt);
            long long nw = old + cnt;
            if(itaw != mp[ru].end()) itaw->second = (int)nw;
            else mp[ru].emplace(w, (int)nw);
            tot2 += C2(nw);

            auto itwrv = mp[w].find(rv);
            if(itwrv != mp[w].end()) mp[w].erase(itwrv);
            auto itwru = mp[w].find(ru);
            if(itwru != mp[w].end()) itwru->second = (int)nw;
            else mp[w].emplace(ru, (int)nw);
        }
        mp[rv].clear();

        tot1 += C2((long long)sz[ru]-f[ru]);
    }

    long long getAns() const { return tot1 - tot2; }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n,m;
    if(!(cin>>n>>m)) return 0;
    DSU dsu(2*n);
    for(int i=0;i<n;i++){
        dsu.add_split_pair(i, n+i);
    }
    vector<int>a(m),b(m);
    for(int i=0;i<m;i++){cin>>a[i]>>b[i]; --a[i]; --b[i];}
    for(int i=0;i<m;i++){
        dsu.unite(a[i], n+b[i]);
        cout<<dsu.getAns()<<"\n";
    }
    return 0;
}