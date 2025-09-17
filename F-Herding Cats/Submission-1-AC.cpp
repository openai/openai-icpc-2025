#include <bits/stdc++.h>
using namespace std;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t;
    if(!(cin >> t)) return 0;
    while(t--){
        int n, m;
        cin >> n >> m;

        vector<int> LB(m+1, 0);
        vector<vector<int>> groups(m+1);
        vector<vector<int>> likes(n);
        vector<int> p(n);

        for(int i=0;i<n;i++){
            int pi, k;
            cin >> pi >> k;
            p[i] = pi;
            groups[pi].push_back(i);
            likes[i].resize(k);
            for(int j=0;j<k;j++){
                int v; cin >> v;
                likes[i][j] = v;
                if (LB[v] < pi) LB[v] = pi;
            }
        }

        // Prefix check: enough plants with LB <= each position
        vector<int> freq(m+1, 0);
        for(int v=1; v<=m; ++v) freq[LB[v]]++;
        long long acc = freq[0];
        bool ok = true;
        for(int j=1; j<=m; ++j){
            acc += freq[j];
            if(acc < j){
                ok = false; break;
            }
        }

        // Check for each position j with cats that intersection contains a plant v with LB[v] == j
        if(ok){
            vector<int> cnt(m+1, 0);
            vector<int> touched;
            for(int j=1; j<=m && ok; ++j){
                auto &g = groups[j];
                int need = (int)g.size();
                if(need == 0) continue;
                touched.clear();
                for(int idx : g){
                    for(int v : likes[idx]){
                        if(cnt[v]==0) touched.push_back(v);
                        cnt[v]++;
                    }
                }
                bool found=false;
                for(int v : touched){
                    if(cnt[v]==need && LB[v]==j){
                        found=true; break;
                    }
                }
                for(int v : touched) cnt[v]=0;
                if(!found) ok=false;
            }
        }

        cout << (ok ? "yes" : "no") << '\n';
    }
    return 0;
}