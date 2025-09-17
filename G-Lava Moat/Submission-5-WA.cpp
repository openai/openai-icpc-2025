#include <bits/stdc++.h>
using namespace std;

struct DSU {
    struct Change{
        int type; //0 same component (cycle), 1 union
        int r1;
        int r2;
        int size_r1;
        double a1,b1,a2,b2;
    };
    vector<int> parent;
    vector<int> sz;
    vector<double> sumA, sumB;
    vector<Change> hist;
    DSU(int n=0){ init(n); }
    void init(int n){
        parent.resize(n);
        sz.assign(n,1);
        sumA.assign(n,0.0);
        sumB.assign(n,0.0);
        for(int i=0;i<n;i++) parent[i]=i;
        hist.clear();
    }
    int find(int x){
        while(parent[x]!=x) x=parent[x];
        return x;
    }
    int snapshot() { return (int)hist.size(); }
    void rollback(int snap){
        while((int)hist.size()>snap){
            Change ch=hist.back(); hist.pop_back();
            if(ch.type==0){
                int r=ch.r1;
                sumA[r]=ch.a1; sumB[r]=ch.b1;
            } else {
                int r=ch.r1, c=ch.r2;
                parent[c]=c;
                sz[r]=ch.size_r1;
                sumA[r]=ch.a1; sumB[r]=ch.b1;
                sumA[c]=ch.a2; sumB[c]=ch.b2;
            }
        }
    }
    void unite(int u,int v,double A,double B){
        int ru=find(u), rv=find(v);
        if(ru==rv){
            Change ch;
            ch.type=0; ch.r1=ru; ch.r2=-1; ch.size_r1=0;
            ch.a1=sumA[ru]; ch.b1=sumB[ru]; ch.a2=ch.b2=0;
            hist.push_back(ch);
            sumA[ru]+=A; sumB[ru]+=B;
        } else {
            if(sz[ru]<sz[rv]) swap(ru,rv);
            Change ch;
            ch.type=1; ch.r1=ru; ch.r2=rv; ch.size_r1=sz[ru];
            ch.a1=sumA[ru]; ch.b1=sumB[ru]; ch.a2=sumA[rv]; ch.b2=sumB[rv];
            hist.push_back(ch);
            parent[rv]=ru;
            sz[ru]+=sz[rv];
            sumA[ru]+=sumA[rv]+A;
            sumB[ru]+=sumB[rv]+B;
        }
    }
};

struct Vertex{
    double x,y,z;
};
struct ContourEdge{
    int u,v;
    int start,end;
    double A,B;
};

static inline long long makeKey(int u,int v,int n){
    if(u>v) swap(u,v);
    return (long long)u*((long long)n+1) + (long long)v;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    if(!(cin>>t)) return 0;
    while(t--){
        int w,l;
        int n,m;
        cin>>w>>l>>n>>m;
        vector<Vertex> verts(n+1);
        int sw=-1,se=-1,nw=-1,ne=-1;
        for(int i=1;i<=n;i++){
            int x,y,z;
            cin>>x>>y>>z;
            verts[i]={double(x),double(y),double(z)};
            if(x==0 && y==0) sw=i;
            else if(x==w && y==0) se=i;
            else if(x==0 && y==l) nw=i;
            else if(x==w && y==l) ne=i;
        }

        unordered_map<long long,int> edgeMap;
        edgeMap.reserve(m*6);
        vector<pair<int,int>> edges;
        edges.reserve(m*3);

        auto getEdgeId = [&](int u,int v)->int{
            long long key=makeKey(u,v,n);
            auto it=edgeMap.find(key);
            if(it==edgeMap.end()){
                int id= (int)edges.size();
                edgeMap[key]=id;
                edges.push_back({u,v});
                return id;
            }
            return it->second;
        };

        vector<array<int,3>> triangles;
        triangles.reserve(m);
        for(int i=0;i<m;i++){
            int a,b,c;
            cin>>a>>b>>c;
            triangles.push_back({a,b,c});
            getEdgeId(a,b); getEdgeId(b,c); getEdgeId(c,a);
        }

        int E=edges.size();

        vector<int> ord(n);
        for(int i=0;i<n;i++) ord[i]=i+1;
        sort(ord.begin(),ord.end(),[&](int a,int b){return verts[a].z<verts[b].z;});
        vector<int> rank(n+1);
        vector<double> sorted_z(n);
        for(int i=0;i<n;i++){
            rank[ord[i]]=i;
            sorted_z[i]=verts[ord[i]].z;
        }

        int westId=-1,eastId=-1;
        {
            long long key=makeKey(sw,nw,n);
            auto it=edgeMap.find(key);
            if(it!=edgeMap.end()) westId=it->second;
            key=makeKey(se,ne,n);
            it=edgeMap.find(key);
            if(it!=edgeMap.end()) eastId=it->second;
        }

        int NInt=n-1;

        vector<ContourEdge> cedges;
        cedges.reserve(m*2);

        for(auto &tri: triangles){
            int ids[3]={tri[0],tri[1],tri[2]};
            sort(ids,ids+3,[&](int a,int b){return verts[a].z<verts[b].z;});
            int low=ids[0], mid=ids[1], high=ids[2];
            int rl=rank[low], rm=rank[mid], rh=rank[high];
            double zl=verts[low].z, zm=verts[mid].z, zh=verts[high].z;
            double xl=verts[low].x, yl=verts[low].y;
            double xm=verts[mid].x, ym=verts[mid].y;
            double xh=verts[high].x, yh=verts[high].y;

            double dh=zh-zl;
            double dm=zm-zl;
            double dmh=zh-zm;

            double dxLow=(xh-xl)/dh - (xm-xl)/dm;
            double dyLow=(yh-yl)/dh - (ym-yl)/dm;
            double Klow=std::sqrt(dxLow*dxLow+dyLow*dyLow);

            double dxHigh=((xm-xh)/dmh) - ((xl-xh)/dh);
            double dyHigh=((ym-yh)/dmh) - ((yl-yh)/dh);
            double Khigh=std::sqrt(dxHigh*dxHigh+dyHigh*dyHigh);

            int idLM=getEdgeId(low,mid);
            int idLH=getEdgeId(low,high);
            int idMH=getEdgeId(mid,high);

            int start=rl;
            int end=rm-1;
            if(start<=end && start<=NInt-1 && end>=0){
                if(end> NInt-1) end=NInt-1;
                ContourEdge e;
                e.u=idLM; e.v=idLH; e.start=start; e.end=end;
                e.A=Klow; e.B=-Klow*zl;
                cedges.push_back(e);
            }

            start=rm; end=rh-1;
            if(start<=end && start<=NInt-1 && end>=0){
                if(end> NInt-1) end=NInt-1;
                ContourEdge e;
                e.u=idMH; e.v=idLH; e.start=start; e.end=end;
                e.A=-Khigh; e.B=Khigh*zh;
                cedges.push_back(e);
            }
        }

        vector<vector<int>> seg;
        if(NInt>0) seg.assign(4*NInt, {});

        function<void(int,int,int,int,int,int)> addSeg = [&](int node,int l_int,int r_int,int ql,int qr,int idx){
            if(ql>r_int||qr<l_int) return;
            if(ql<=l_int && r_int<=qr){
                seg[node].push_back(idx);
                return;
            }
            int mid=(l_int+r_int)/2;
            addSeg(node*2,l_int,mid,ql,qr,idx);
            addSeg(node*2+1,mid+1,r_int,ql,qr,idx);
        };

        for(int i=0;i<(int)cedges.size();i++){
            if(NInt<=0) continue;
            int s=cedges[i].start, eidx=cedges[i].end;
            if(s>eidx) continue;
            if(s<0) s=0;
            if(eidx>=NInt) eidx=NInt-1;
            if(s>eidx) continue;
            addSeg(1,0,NInt-1,s,eidx,i);
        }

        DSU dsu(E);
        double ans=1e100;
        bool possible=false;

        function<void(int,int,int)> dfs = [&](int node,int l_int,int r_int){
            int snap=dsu.snapshot();
            for(int idx: seg[node]){
                auto &ce=cedges[idx];
                dsu.unite(ce.u,ce.v,ce.A,ce.B);
            }
            if(l_int==r_int){
                int rw=dsu.find(westId), re=dsu.find(eastId);
                if(rw==re){
                    double A=dsu.sumA[rw], B=dsu.sumB[rw];
                    double cLow=sorted_z[l_int];
                    double cHigh=sorted_z[l_int+1];
                    double valLow=A*cLow+B;
                    double valHigh=A*cHigh+B;
                    double cand=min(valLow,valHigh);
                    if(cand<ans){ ans=cand; possible=true; }
                }
            } else {
                int mid=(l_int+r_int)/2;
                dfs(node*2,l_int,mid);
                dfs(node*2+1,mid+1,r_int);
            }
            dsu.rollback(snap);
        };

        if(NInt>0) dfs(1,0,NInt-1);

        if(!possible){
            cout<<"impossible";
        } else {
            cout.setf(std::ios::fixed);
            cout<<setprecision(9)<<ans;
        }
        if(t) cout<<"\n";
    }
    return 0;
}