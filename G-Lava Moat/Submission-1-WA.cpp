#include <bits/stdc++.h>
using namespace std;

struct DSURollback {
    int n;
    vector<int> parent, sz;
    vector<double> A, B; // total coefficients for edges in component
    struct Change {int a,b; int size_a; double A_a,B_a; bool merged;};
    vector<Change> st;
    DSURollback(int n=0){init(n);}    
    void init(int n_){n=n_; parent.resize(n); sz.assign(n,1); A.assign(n,0.0); B.assign(n,0.0); st.clear(); for(int i=0;i<n;i++) parent[i]=i;}
    int find(int x){ while(parent[x]!=x) x=parent[x]; return x; }
    void unite(int u,int v,double a,double b){
        int ru=find(u), rv=find(v);
        if(ru==rv){ st.push_back({0,0,0,0,0,false}); return; }
        if(sz[ru]<sz[rv]) swap(ru,rv);
        st.push_back({ru,rv,sz[ru],A[ru],B[ru],true});
        parent[rv]=ru; sz[ru]+=sz[rv];
        A[ru]+=A[rv]+a; B[ru]+=B[rv]+b; 
    }
    int snapshot(){ return (int)st.size(); }
    void rollback(int snap){
        while((int)st.size()>snap){
            Change c=st.back(); st.pop_back();
            if(!c.merged) continue; 
            int ru=c.a, rv=c.b; parent[rv]=rv; sz[ru]=c.size_a; A[ru]=c.A_a; B[ru]=c.B_a; 
        }
    }
};

struct Conn{int u,v; double A,B;};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t; 
    if(!(cin>>t)) return 0;
    cout.setf(std::ios::fixed);
    cout<<setprecision(10);
    for(int tc=0; tc<t; ++tc){
        long long W,L; int n,m; 
        cin>>W>>L>>n>>m; 
        vector<double> x(n+1), y(n+1), z(n+1);
        for(int i=1;i<=n;i++){ long long xi,yi,zi; cin>>xi>>yi>>zi; x[i]=xi; y[i]=yi; z[i]=zi; }
        vector<array<int,3>> tri(m);
        for(int i=0;i<m;i++){ int a,b,c; cin>>a>>b>>c; tri[i]={a,b,c}; }
        // Build edges map
        unordered_map<long long,int> mp; mp.reserve(m*3*2); mp.max_load_factor(0.7);
        auto keypair=[&](int a,int b){ if(a>b) swap(a,b); return ((long long)a<<32) ^ (unsigned long long)b; };
        struct Edge{int u,v;};
        vector<Edge> edges; edges.reserve(3*m);
        auto getEdge=[&](int a,int b){ 
            long long k=keypair(a,b); auto it=mp.find(k); if(it==mp.end()){int id=edges.size(); edges.push_back({min(a,b),max(a,b)}); mp[k]=id; return id;} else return it->second; };
        for(int i=0;i<m;i++){
            int a=tri[i][0], b=tri[i][1], c=tri[i][2];
            getEdge(a,b); getEdge(b,c); getEdge(c,a);
        }
        int E=edges.size();
        
        // identify corners: (0,0), (W,0), (W,L), (0,L)
        int idSW=-1,idSE=-1,idNE=-1,idNW=-1; 
        for(int i=1;i<=n;i++){
            if(x[i]==0 && y[i]==0) idSW=i; 
            if(x[i]==W && y[i]==0) idSE=i; 
            if(x[i]==W && y[i]==L) idNE=i; 
            if(x[i]==0 && y[i]==L) idNW=i; 
        }
        int westEdge = getEdge(idSW,idNW);
        int eastEdge = getEdge(idSE,idNE);
        
        // Sort heights
        vector<double> H; H.reserve(n); for(int i=1;i<=n;i++) H.push_back(z[i]); sort(H.begin(),H.end()); 
        int K=H.size();
        vector<int> pos(n+1);
        map<long long,int> mapInt; for(int i=0;i<K;i++){ mapInt[(long long)H[i]]=i; }
        for(int i=1;i<=n;i++){ pos[i] = mapInt[(long long)z[i]]; }
        
        int numTimes = K-1; 
        if(numTimes<=0){
            cout << "impossible\n";
            continue;
        }
        vector<vector<Conn>> seg(4*max(1,numTimes));
        function<void(int,int,int,int,int,const Conn&)> addInterval = [&](int node, int l, int r, int Lq, int Rq, const Conn& e){
            if(Lq>=Rq || l>=Rq || r<=Lq) return; 
            if(Lq<=l && r<=Rq){ seg[node].push_back(e); return; }
            int mid=(l+r)/2; 
            addInterval(node*2,l,mid,Lq,Rq,e); 
            addInterval(node*2+1,mid,r,Lq,Rq,e);
        };
        
        auto edgeId=[&](int a,int b){ return getEdge(a,b); };
        auto makeS = [&](int u,int v,int w){
            double dzuv = z[v]-z[u]; double dzuw = z[w]-z[u];
            pair<double,double> v1 = { (x[v]-x[u])/dzuv, (y[v]-y[u])/dzuv };
            pair<double,double> v2 = { (x[w]-x[u])/dzuw, (y[w]-y[u])/dzuw };
            double dx = v1.first - v2.first; double dy = v1.second - v2.second; return sqrt(dx*dx + dy*dy);
        };
        // Build connection intervals for each triangle
        for(int i=0;i<m;i++){
            int a=tri[i][0], b=tri[i][1], c=tri[i][2];
            int i0=a, i1=b, i2=c; 
            if(z[i1] < z[i0]) swap(i0,i1); if(z[i2] < z[i1]) swap(i1,i2); if(z[i1] < z[i0]) swap(i0,i1); 
            int e1 = edgeId(i0,i1), e2 = edgeId(i0,i2);
            double S1 = makeS(i0,i1,i2);
            int L1 = pos[i0], R1 = pos[i1];
            if(L1 < R1){ Conn e{e1,e2,+S1,-S1*z[i0]}; addInterval(1,0,numTimes,L1,R1,e); }
            int e3 = edgeId(i2,i0), e4 = edgeId(i2,i1);
            double S2 = makeS(i2,i0,i1);
            int L2 = pos[i1], R2 = pos[i2];
            if(L2 < R2){ Conn e{e3,e4,-S2, S2*z[i2]}; addInterval(1,0,numTimes,L2,R2,e); }
        }
        
        // intersection of west and east active intervals
        int Lw = min(pos[idSW], pos[idNW]);
        int Rw = max(pos[idSW], pos[idNW]);
        int Le = min(pos[idSE], pos[idNE]);
        int Re = max(pos[idSE], pos[idNE]);
        int Lint = max(Lw, Le); int Rint = min(Rw, Re); 
        if(Rint<=Lint){
            cout << "impossible\n";
            continue;
        }
        double startBound = H[Lint];
        double endBound = H[Rint];
        
        DSURollback dsu(E);
        double best = 1e300; bool any=false;
        function<void(int,int,int)> dfs = [&](int node, int l, int r){
            int snap=dsu.snapshot();
            for(const Conn &e: seg[node]) dsu.unite(e.u, e.v, e.A, e.B);
            if(l+1==r){
                int j=l; double left=H[j], right=H[j+1];
                double Lb = max(left, startBound); double Rb = min(right, endBound);
                if(Lb < Rb){
                    int rw = dsu.find(westEdge), re = dsu.find(eastEdge);
                    if(rw==re){
                        double A = dsu.A[rw], B = dsu.B[rw];
                        double v1 = A*Lb + B; double v2 = A*Rb + B; double val=min(v1,v2); 
                        if(val<best){best=val; any=true;}
                    }
                }
            } else {
                int mid=(l+r)/2; dfs(node*2,l,mid); dfs(node*2+1,mid,r); 
            }
            dsu.rollback(snap);
        };
        dfs(1,0,numTimes);
        if(!any) cout << "impossible\n"; else cout << best << "\n";
    }
    return 0;
}