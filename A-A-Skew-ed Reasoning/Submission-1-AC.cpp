#include <bits/stdc++.h>
using namespace std;

struct Task {
    int u; bool isMin; long long base, step; long long start, len;
};

long long ceil_div_ll(long long a, long long b){
    if(a>=0) return (a + b - 1)/b;
    else return - ((-a)/b);
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int n;
    if(!(cin>>n)) return 0;
    vector<int> L(n+1), R(n+1);
    for(int i=1;i<=n;i++){
        cin>>L[i]>>R[i];
    }

    // compute subtree sizes using topological order (children indices > parent)
    vector<int> sz(n+1,0);
    for(int i=n;i>=1;i--){
        sz[i] = 1 + (L[i]?sz[L[i]]:0) + (R[i]?sz[R[i]]:0);
    }

    vector<char> minA(n+1), maxA(n+1); // whether A (stream starting empty) is left child
    bool ok=true;
    for(int u=1; u<=n; ++u){
        int nl = L[u]?sz[L[u]]:0;
        int nr = R[u]?sz[R[u]]:0;
        if(nl==0 && nr>0){ ok=false; break; }
        bool opt1 = (nl>=1 && nl<=nr+1); // A=Left feasible
        bool opt2 = (nl>=nr);            // A=Right feasible
        if(!opt1 && !opt2){ ok=false; break; }
        if(opt1 && opt2){
            if(nl==nr+1){ minA[u]=1; maxA[u]=0; }
            else if(nl==nr){ minA[u]=0; maxA[u]=1; }
            else { minA[u]=0; maxA[u]=1; } // should not occur
        } else if(opt1){
            minA[u]=1; maxA[u]=1;
        } else { // opt2
            minA[u]=0; maxA[u]=0;
        }
    }
    if(!ok){
        cout<<"impossible\n";
        return 0;
    }

    auto fill_mode = [&](bool isMin){
        vector<int> res(n,0);
        vector<Task> st; st.reserve(2*n+10);
        st.push_back({1,isMin,0,1,0,(long long)sz[1]});
        while(!st.empty()){
            Task t = st.back(); st.pop_back();
            int u=t.u; if(u==0 || t.len<=0) continue;
            long long base=t.base, step=t.step, start=t.start, len=t.len, end=start+len;
            int l=L[u], r=R[u];
            if(l==0 && r==0){
                res[base]=u;
                continue;
            }
            bool ALeft = isMin ? minA[u] : maxA[u];
            int A = ALeft? l : r;
            int B = ALeft? r : l;
            int nA = A? sz[A] : 0;
            int nB = B? sz[B] : 0;
            long long k = ALeft ? (long long)(nB - nA + 1) : (long long)(nB - nA);

            // B prefix
            long long bL = max(0LL, start), bR = min(end, k);
            if(B && bL < bR){
                st.push_back({B, isMin, base + step*(bL - start), step, bL, bR-bL});
            }
            // root
            if(start <= k && k < end){
                res[base + step*(k - start)] = u;
            }
            // A zipped
            long long lowA = max(0LL, ceil_div_ll(start - (k+1), 2));
            long long highA = min((long long)nA, ceil_div_ll(end - (k+1), 2));
            if(A && lowA < highA){
                st.push_back({A, isMin, base + step*((k+1 - start) + 2*lowA), step*2, lowA, highA-lowA});
            }
            // B remainder zipped
            long long remB = nB - k; if(remB < 0) remB = 0;
            long long lowB = max(0LL, ceil_div_ll(start - (k+2), 2));
            long long highB = min(remB, ceil_div_ll(end - (k+2), 2));
            if(B && lowB < highB){
                st.push_back({B, isMin, base + step*((k+2 - start) + 2*lowB), step*2, k + lowB, highB - lowB});
            }
        }
        return res;
    };

    vector<int> resMin = fill_mode(true);
    vector<int> resMax = fill_mode(false);

    for(int i=0;i<n;i++){ if(i) cout<<' '; cout<<resMin[i]; } cout<<"\n";
    for(int i=0;i<n;i++){ if(i) cout<<' '; cout<<resMax[i]; } cout<<"\n";
    return 0;
}