#include <bits/stdc++.h>
using namespace std;

struct DSU {
    int n;
    vector<int> parent, sz;
    vector<double> A, B; // aggregated coefficients per component
    struct Change {
        int type; // 1 union, 2 add within component
        int a, b; // for union: roots; for add: root
        int prev_sz; // size of a before union
        double oldA_a, oldB_a; // previous A,B of root a before union
        double addA, addB; // edge coefficients added
    };
    vector<Change> st;

    DSU(int n=0){ init(n); }
    void init(int n_){
        n = n_;
        parent.resize(n);
        sz.assign(n,1);
        A.assign(n,0.0);
        B.assign(n,0.0);
        for(int i=0;i<n;i++) parent[i]=i;
        st.clear();
    }
    int find(int x){
        while(parent[x]!=x) x=parent[x];
        return x;
    }
    void add_edge(int u,int v,double a,double b){
        int ru=find(u), rv=find(v);
        if(ru==rv){
            st.push_back({2,ru,-1,0,0,0,a,b});
            A[ru]+=a; B[ru]+=b;
            return;
        }
        if(sz[ru]<sz[rv]) swap(ru,rv);
        st.push_back({1,ru,rv,sz[ru],A[ru],B[ru],a,b});
        parent[rv]=ru; sz[ru]+=sz[rv];
        A[ru]=A[ru]+A[rv]+a; B[ru]=B[ru]+B[rv]+b;
    }
    void rollback(int target){
        while((int)st.size()>target){
            Change ch=st.back(); st.pop_back();
            if(ch.type==2){
                A[ch.a]-=ch.addA; B[ch.a]-=ch.addB;
            }else{
                int ru=ch.a, rv=ch.b;
                parent[rv]=rv;
                sz[ru]=ch.prev_sz;
                A[ru]=ch.oldA_a; B[ru]=ch.oldB_a;
            }
        }
    }
};

struct Event{ int u,v; double a,b; };

struct SegTree{
    int n;
    vector<vector<Event>> tree;
    void init(int n_){ n=n_; tree.assign(4*n+5,{}); }
    void add(int node,int l,int r,int ql,int qr, const Event &e){
        if(ql>qr) return;
        if(ql<=l && r<=qr){ tree[node].push_back(e); return; }
        int mid=(l+r)/2;
        if(ql<=mid) add(node*2,l,mid,ql,qr,e);
        if(qr>mid) add(node*2+1,mid+1,r,ql,qr,e);
    }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t; 
    if(!(cin>>t)) return 0;
    cout.setf(std::ios::fixed); cout<<setprecision(10);
    while(t--){
        long long w,l; int n,m;
        cin>>w>>l>>n>>m;
        vector<long long> X(n+1), Y(n+1), Z(n+1);
        for(int i=1;i<=n;i++) cin>>X[i]>>Y[i]>>Z[i];
        vector<array<int,3>> tris(m);
        for(int i=0;i<m;i++){ int a,b,c; cin>>a>>b>>c; tris[i]={a,b,c}; }

        int id_sw=-1,id_nw=-1,id_se=-1,id_ne=-1;
        for(int i=1;i<=n;i++){
            if(X[i]==0 && Y[i]==0) id_sw=i;
            else if(X[i]==0 && Y[i]==l) id_nw=i;
            else if(X[i]==w && Y[i]==0) id_se=i;
            else if(X[i]==w && Y[i]==l) id_ne=i;
        }

        struct Key{int a,b;};
        struct Hash{ size_t operator()(Key const&k) const noexcept { return ((uint64_t)k.a<<32) ^ k.b; } };
        struct Eq{ bool operator()(Key const&x, Key const&y) const noexcept { return x.a==y.a && x.b==y.b; } };
        unordered_map<Key,int,Hash,Eq> edgeId;
        edgeId.reserve(m*6+10);
        edgeId.max_load_factor(0.7);
        int E=0;
        auto getEdge=[&](int u,int v){
            if(u>v) swap(u,v);
            Key k{u,v};
            auto it=edgeId.find(k);
            if(it==edgeId.end()){ edgeId[k]=E++; return E-1; }
            else return it->second;
        };
        for(auto &tr: tris){
            int a=tr[0],b=tr[1],c=tr[2];
            getEdge(a,b); getEdge(b,c); getEdge(c,a);
        }
        int westEdgeId=getEdge(id_sw,id_nw);
        int eastEdgeId=getEdge(id_se,id_ne);

        vector<long long> zvals; zvals.reserve(n);
        for(int i=1;i<=n;i++) zvals.push_back(Z[i]);
        sort(zvals.begin(), zvals.end());
        zvals.erase(unique(zvals.begin(), zvals.end()), zvals.end());
        int K=zvals.size();
        if(K<=1){
            cout<<"impossible\n";
            continue;
        }
        auto idx=[&](long long z){ return (int)(lower_bound(zvals.begin(), zvals.end(), z) - zvals.begin()); };
        SegTree seg; seg.init(K-1);
        auto vec2=[&](int u,int v){
            return pair<double,double>{ (double)(X[v]-X[u]), (double)(Y[v]-Y[u]) };
        };
        for(auto &tr: tris){
            int a=tr[0],b=tr[1],c=tr[2];
            int ids[3]={a,b,c};
            long long Zs[3]={Z[a],Z[b],Z[c]};
            int ord[3]={0,1,2};
            sort(ord,ord+3,[&](int i,int j){ return Zs[i] < Zs[j]; });
            int i1=ids[ord[0]], i2=ids[ord[1]], i3=ids[ord[2]];
            long long z1=Z[i1], z2=Z[i2], z3=Z[i3];
            int e12=getEdge(i1,i2), e13=getEdge(i1,i3), e23=getEdge(i2,i3);
            auto v12=vec2(i1,i2); auto v13=vec2(i1,i3);
            double inv12 = 1.0/(double)(z2 - z1);
            double inv13 = 1.0/(double)(z3 - z1);
            double dx = v12.first*inv12 - v13.first*inv13;
            double dy = v12.second*inv12 - v13.second*inv13;
            double c_low = sqrt(dx*dx + dy*dy);
            auto v31=vec2(i3,i1); auto v32=vec2(i3,i2);
            double inv31 = 1.0/(double)(z1 - z3);
            double inv32 = 1.0/(double)(z2 - z3);
            dx = v31.first*inv31 - v32.first*inv32;
            dy = v31.second*inv31 - v32.second*inv32;
            double c_high = sqrt(dx*dx + dy*dy);

            int L1=idx(z1), R1=idx(z2)-1;
            if(L1<=R1){
                Event e{e12, e13, c_low, -c_low*(double)z1};
                seg.add(1,0,K-2,L1,R1,e);
            }
            int L2=idx(z2), R2=idx(z3)-1;
            if(L2<=R2){
                Event e{e13, e23, -c_high, c_high*(double)z3};
                seg.add(1,0,K-2,L2,R2,e);
            }
        }

        vector<char> activeW(K-1,0), activeE(K-1,0);
        int wLidx = idx(min(Z[id_sw], Z[id_nw]));
        int wRidx = idx(max(Z[id_sw], Z[id_nw])) - 1;
        for(int j=max(0,wLidx); j<=min(K-2,wRidx); ++j) activeW[j]=1;
        int eLidx = idx(min(Z[id_se], Z[id_ne]));
        int eRidx = idx(max(Z[id_se], Z[id_ne])) - 1;
        for(int j=max(0,eLidx); j<=min(K-2,eRidx); ++j) activeE[j]=1;

        DSU dsu(E);
        double best = 1e100;
        bool found = false;

        function<void(int,int,int)> dfs = [&](int node,int l,int r){
            int save = dsu.st.size();
            for(const Event &ev : seg.tree[node]){
                dsu.add_edge(ev.u, ev.v, ev.a, ev.b);
            }
            if(l==r){
                if(activeW[l] && activeE[l]){
                    int rw = dsu.find(westEdgeId);
                    int re = dsu.find(eastEdgeId);
                    if(rw == re){
                        double A = dsu.A[rw], B = dsu.B[rw];
                        double left = A * (double)zvals[l] + B;
                        double right = A * (double)zvals[l+1] + B;
                        double val = min(left, right);
                        if(isfinite(val)){
                            best = min(best, val);
                            found = true;
                        }
                    }
                }
                dsu.rollback(save);
                return;
            }
            int mid=(l+r)/2;
            dfs(node*2,l,mid);
            dfs(node*2+1,mid+1,r);
            dsu.rollback(save);
        };

        if(K-1>0) dfs(1,0,K-2);
        if(!found) cout<<"impossible\n";
        else cout<<best<<"\n";
    }
    return 0;
}