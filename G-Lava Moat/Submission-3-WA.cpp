#include <bits/stdc++.h>
using namespace std;

struct IntervalEdge{
    int u,v;
    int l,r; // inclusive gap indices
    double a,b;
};

struct DSU{
    vector<int> parent, sz;
    vector<double> suma, sumb;
    struct Change{
        int type; // 0 union, 1 internal add
        int x;
        int prev_sz;
        double prev_a;
        double prev_b;
    };
    vector<Change> st;
    DSU(int n=0){ init(n); }
    void init(int n){
        parent.resize(n);
        sz.assign(n,1);
        suma.assign(n,0);
        sumb.assign(n,0);
        st.clear();
        for(int i=0;i<n;i++) parent[i]=i;
    }
    int find(int x){
        while(parent[x]!=x) x=parent[x];
        return x;
    }
    int snap() const { return (int)st.size(); }
    void rollback(int snap){
        while((int)st.size()>snap){
            auto ch=st.back(); st.pop_back();
            if(ch.type==0){
                int rv=ch.x;
                int ru=parent[rv];
                parent[rv]=rv;
                sz[ru]=ch.prev_sz;
                suma[ru]=ch.prev_a;
                sumb[ru]=ch.prev_b;
            }else{
                int r=ch.x;
                suma[r]=ch.prev_a;
                sumb[r]=ch.prev_b;
            }
        }
    }
    void unite(int u,int v,double a,double b){
        int ru=find(u), rv=find(v);
        if(ru!=rv){
            if(sz[ru]<sz[rv]) swap(ru,rv);
            st.push_back({0, rv, sz[ru], suma[ru], sumb[ru]});
            parent[rv]=ru;
            sz[ru]+=sz[rv];
            suma[ru]+=suma[rv]+a;
            sumb[ru]+=sumb[rv]+b;
        }else{
            st.push_back({1, ru, 0, suma[ru], sumb[ru]});
            suma[ru]+=a;
            sumb[ru]+=b;
        }
    }
};

struct Tri{int a,b,c;};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    if(!(cin>>t)) return 0;
    while(t--){
        int w,l;
        int n,m;
        cin>>w>>l>>n>>m;
        vector<double> vx(n+1), vy(n+1), vz(n+1);
        int SW=-1, SE=-1, NW=-1, NE=-1;
        for(int i=1;i<=n;i++){
            int xi, yi, zi;
            cin>>xi>>yi>>zi;
            vx[i]=xi; vy[i]=yi; vz[i]=zi;
            if(xi==0 && yi==0) SW=i;
            else if(xi==w && yi==0) SE=i;
            else if(xi==w && yi==l) NE=i;
            else if(xi==0 && yi==l) NW=i;
        }
        vector<Tri> tris(m);
        for(int i=0;i<m;i++){
            int a,b,c;
            cin>>a>>b>>c;
            tris[i]={a,b,c};
        }

        vector<int> ord(n);
        iota(ord.begin(),ord.end(),1);
        sort(ord.begin(),ord.end(),[&](int i,int j){return vz[i]<vz[j];});
        vector<int> rank(n+1);
        vector<double> Z(n+1);
        for(int idx=0;idx<n;idx++){
            rank[ord[idx]]=idx+1;
            Z[idx+1]=vz[ord[idx]];
        }

        double westMin=min(vz[SW],vz[NW]), westMax=max(vz[SW],vz[NW]);
        double eastMin=min(vz[SE],vz[NE]), eastMax=max(vz[SE],vz[NE]);
        double I_low=max(westMin,eastMin);
        double I_high=min(westMax,eastMax);

        if(!(I_low < I_high)){
            cout<<"impossible\n";
            continue;
        }

        long long base = (long long)n + 1;
        unordered_map<long long,int> edgeMap;
        edgeMap.reserve(m*6+10);
        vector<pair<int,int>> edges;
        auto getEdge=[&](int u,int v)->int{
            if(u>v) swap(u,v);
            long long key = (long long)u*base + v;
            auto it=edgeMap.find(key);
            if(it!=edgeMap.end()) return it->second;
            int id=(int)edges.size();
            edges.push_back({u,v});
            edgeMap[key]=id;
            return id;
        };

        for(auto &tr: tris){
            getEdge(tr.a,tr.b);
            getEdge(tr.b,tr.c);
            getEdge(tr.c,tr.a);
        }

        int E=edges.size();
        int westEdge=getEdge(SW,NW);
        int eastEdge=getEdge(SE,NE);

        vector<IntervalEdge> intervalEdges;
        intervalEdges.reserve(m*2+10);

        for(auto &tr: tris){
            int ids[3]={tr.a,tr.b,tr.c};
            int p[3]={0,1,2};
            sort(p,p+3,[&](int ii,int jj){return vz[ids[ii]]<vz[ids[jj]];});
            int A=ids[p[0]], B=ids[p[1]], C=ids[p[2]];
            double zA=vz[A], zB=vz[B], zC=vz[C];
            int rA=rank[A], rB=rank[B], rC=rank[C];
            double xA=vx[A], yA=vy[A];
            double xB=vx[B], yB=vy[B];
            double xC=vx[C], yC=vy[C];
            int eAB=getEdge(A,B);
            int eAC=getEdge(A,C);
            int eBC=getEdge(B,C);

            int lIdx=rA;
            int rIdx=rB-1;
            if(lIdx<=rIdx){
                double invCA=1.0/(zC - zA);
                double invBA=1.0/(zB - zA);
                double dx = (xC - xA)*invCA - (xB - xA)*invBA;
                double dy = (yC - yA)*invCA - (yB - yA)*invBA;
                double K_low = hypot(dx,dy);
                double a = K_low;
                double b = -zA*K_low;
                intervalEdges.push_back({eAB,eAC,lIdx,rIdx,a,b});
            }

            lIdx=rB;
            rIdx=rC-1;
            if(lIdx<=rIdx){
                double invCA=1.0/(zC - zA);
                double invCB=1.0/(zC - zB);
                double ux = (xC - xA)*invCA;
                double uy = (yC - yA)*invCA;
                double wx = (xC - xB)*invCB;
                double wy = (yC - yB)*invCB;
                double dx = ux - wx;
                double dy = uy - wy;
                double K_high = hypot(dx,dy);
                double a = -K_high;
                double b = zC*K_high;
                intervalEdges.push_back({eAC,eBC,lIdx,rIdx,a,b});
            }
        }

        int T = n-1;
        vector<vector<int>> seg(4*T+10);

        function<void(int,int,int,int,int,int)> addSeg = [&](int idx,int l,int r,int ql,int qr,int id){
            if(ql>r||qr<l) return;
            if(ql<=l&&r<=qr){
                seg[idx].push_back(id);
                return;
            }
            int mid=(l+r)/2;
            addSeg(idx*2,l,mid,ql,qr,id);
            addSeg(idx*2+1,mid+1,r,ql,qr,id);
        };

        for(int i=0;i<(int)intervalEdges.size();i++){
            auto &ie=intervalEdges[i];
            int ll=ie.l, rr=ie.r;
            if(ll>rr) continue;
            if(ll<1) ll=1;
            if(rr>T) rr=T;
            if(ll<=rr) addSeg(1,1,T,ll,rr,i);
        }

        DSU dsu(E);
        double best=numeric_limits<double>::infinity();

        function<void(int,int,int)> dfs = [&](int idx,int l,int r){
            int snap=dsu.snap();
            for(int eid: seg[idx]){
                auto &ie=intervalEdges[eid];
                dsu.unite(ie.u, ie.v, ie.a, ie.b);
            }
            if(l==r){
                double start=Z[l];
                double end=Z[l+1];
                double overlapL=max(start, I_low);
                double overlapR=min(end, I_high);
                if(overlapL<=overlapR){
                    int rw=dsu.find(westEdge);
                    int re=dsu.find(eastEdge);
                    if(rw==re){
                        double A=dsu.suma[rw];
                        double B=dsu.sumb[rw];
                        double vL=A*overlapL+B;
                        double vR=A*overlapR+B;
                        double v=min(vL,vR);
                        if(v<best) best=v;
                    }
                }
            }else{
                int mid=(l+r)/2;
                dfs(idx*2,l,mid);
                dfs(idx*2+1,mid+1,r);
            }
            dsu.rollback(snap);
        };

        if(T>0) dfs(1,1,T);

        if(std::isinf(best)){
            cout<<"impossible\n";
        }else{
            cout.setf(std::ios::fixed);
            cout<<setprecision(9)<<best<<"\n";
        }
    }
    return 0;
}