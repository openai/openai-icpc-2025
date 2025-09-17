#include <bits/stdc++.h>
using namespace std;

struct TempSeg{
    int u,v;
    long long start,end;
    double p,q;
};
struct EdgeItem{
    int u,v;
    double p,q;
};
struct Change{
    int kind;
    int a;
    int b;
    int size_a;
    double p;
    double q;
};
struct RollbackDSU{
    vector<int> parent, sz;
    vector<double> sump, sumq;
    vector<Change> st;
    RollbackDSU(int n=0){ init(n); }
    void init(int n){
        parent.resize(n);
        sz.assign(n,1);
        sump.assign(n,0.0);
        sumq.assign(n,0.0);
        for(int i=0;i<n;i++) parent[i]=i;
        st.clear();
    }
    int find(int x){
        while(parent[x]!=x) x=parent[x];
        return x;
    }
    void addEdge(int u,int v,double p,double q){
        int a=find(u), b=find(v);
        if(a==b){
            sump[a]+=p;
            sumq[a]+=q;
            st.push_back({2,a,-1,0,p,q});
        }else{
            if(sz[a]<sz[b]) swap(a,b);
            st.push_back({1,a,b,sz[a],sump[a],sumq[a]});
            parent[b]=a;
            sz[a]+=sz[b];
            sump[a]+=sump[b]+p;
            sumq[a]+=sumq[b]+q;
        }
    }
    void rollback(int checkpoint){
        while((int)st.size()>checkpoint){
            Change c=st.back(); st.pop_back();
            if(c.kind==2){
                sump[c.a]-=c.p;
                sumq[c.a]-=c.q;
            }else{
                int a=c.a, b=c.b;
                parent[b]=b;
                sz[a]=c.size_a;
                sump[a]=c.p;
                sumq[a]=c.q;
            }
        }
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    if(!(cin>>t)) return 0;
    cout.setf(std::ios::fixed); cout<<setprecision(9);
    while(t--){
        long long w,l;
        int n,m;
        cin>>w>>l>>n>>m;
        vector<long long> x(n), y(n), z(n);
        int sw=-1,se=-1,ne=-1,nw=-1;
        for(int i=0;i<n;i++){
            cin>>x[i]>>y[i]>>z[i];
            if(x[i]==0 && y[i]==0) sw=i;
            else if(x[i]==w && y[i]==0) se=i;
            else if(x[i]==w && y[i]==l) ne=i;
            else if(x[i]==0 && y[i]==l) nw=i;
        }

        vector<array<int,3>> tris;
        tris.reserve(m);

        unordered_map<long long,int> edgeId;
        edgeId.reserve(m*6);
        vector<pair<int,int>> edgeList;
        edgeList.reserve(m*3);
        auto getEdgeId = [&](int a,int b)->int{
            if(a>b) swap(a,b);
            long long key = ((long long)a<<32) | (unsigned int)b;
            auto it=edgeId.find(key);
            if(it!=edgeId.end()) return it->second;
            int id=edgeList.size();
            edgeList.push_back({a,b});
            edgeId[key]=id;
            return id;
        };

        for(int i=0;i<m;i++){
            int a,b,c;
            cin>>a>>b>>c;
            --a;--b;--c;
            tris.push_back({a,b,c});
            getEdgeId(a,b); getEdgeId(b,c); getEdgeId(c,a);
        }

        long long zsw=z[sw], znw=z[nw], zse=z[se], zne=z[ne];
        long long westLow=min(zsw,znw), westHigh=max(zsw,znw);
        long long eastLow=min(zse,zne), eastHigh=max(zse,zne);
        long long overlap_low=max(westLow, eastLow);
        long long overlap_high=min(westHigh, eastHigh);

        if(!(overlap_low < overlap_high)){
            cout<<"impossible\n";
            continue;
        }

        vector<long long> H;
        H.reserve(n+2);
        for(int i=0;i<n;i++) H.push_back(z[i]);
        H.push_back(overlap_low);
        H.push_back(overlap_high);
        sort(H.begin(), H.end());
        H.erase(unique(H.begin(), H.end()), H.end());
        int K=H.size();
        int N=K-1;
        if(N<=0){
            cout<<"impossible\n";
            continue;
        }

        unordered_map<long long,int> idxByH;
        idxByH.reserve(K*2);
        for(int i=0;i<K;i++) idxByH[H[i]]=i;

        int westEdgeId=getEdgeId(sw,nw);
        int eastEdgeId=getEdgeId(se,ne);

        vector<TempSeg> segments;
        segments.reserve(m*2);

        auto addSub=[&](int u,int v,int uu,int vv,long long startH,long long endH,double dx_unit,double dy_unit){
            double dz1=(double)z[v]- (double)z[u];
            double dx1=(double)x[v]- (double)x[u];
            double dy1=(double)y[v]- (double)y[u];
            double ax1=dx1/dz1;
            double ay1=dy1/dz1;
            double bx1=(double)x[u]-ax1*(double)z[u];
            double by1=(double)y[u]-ay1*(double)z[u];

            double dz2=(double)z[vv]- (double)z[uu];
            double dx2=(double)x[vv]- (double)x[uu];
            double dy2=(double)y[vv]- (double)y[uu];
            double ax2=dx2/dz2;
            double ay2=dy2/dz2;
            double bx2=(double)x[uu]-ax2*(double)z[uu];
            double by2=(double)y[uu]-ay2*(double)z[uu];

            double p_raw = dx_unit*(ax2-ax1) + dy_unit*(ay2-ay1);
            double q_raw = dx_unit*(bx2-bx1) + dy_unit*(by2-by1);

            double repZ = ( (double)startH + (double)endH )*0.5;
            double f = p_raw*repZ + q_raw;
            int sign = (f>=0.0 ? 1 : -1);
            if(f==0.0) sign=1;
            double pp=sign*p_raw;
            double qq=sign*q_raw;

            int id1=getEdgeId(u,v);
            int id2=getEdgeId(uu,vv);

            segments.push_back({id1,id2,startH,endH,pp,qq});
        };

        for(auto &tr: tris){
            int a=tr[0], b=tr[1], c=tr[2];
            double x0=x[a], y0=y[a], z0d=z[a];
            double x1=x[b], y1=y[b], z1d=z[b];
            double x2=x[c], y2=y[c], z2d=z[c];
            double ux=x1-x0, uy=y1-y0, uz=z1d-z0d;
            double vx=x2-x0, vy=y2-y0, vz=z2d-z0d;
            double nx= uy*vz - uz*vy;
            double ny= uz*vx - ux*vz;
            double nz= ux*vy - uy*vx;
            double a_grad=-nx/nz;
            double b_grad=-ny/nz;
            double norm = hypot(a_grad,b_grad);
            double dx_unit = -b_grad/norm;
            double dy_unit = a_grad/norm;

            int ids[3]={a,b,c};
            sort(ids, ids+3, [&](int i1,int i2){ return z[i1] < z[i2]; });
            int low=ids[0], mid=ids[1], high=ids[2];
            long long zlow=z[low], zmid=z[mid], zhigh=z[high];

            if(zlow < zmid){
                addSub(low,mid, low,high, zlow,zmid, dx_unit, dy_unit);
            }
            if(zmid < zhigh){
                addSub(low,high, mid,high, zmid,zhigh, dx_unit, dy_unit);
            }
        }

        int sizeTree=1;
        while(sizeTree < N) sizeTree<<=1;
        vector<vector<EdgeItem>> segtree(sizeTree*2);

        for(auto &sg: segments){
            int lidx = idxByH[sg.start];
            int ridx = idxByH[sg.end];
            int L=max(lidx,0);
            int R=min(ridx, N);
            if(L>=R) continue;
            EdgeItem item{sg.u, sg.v, sg.p, sg.q};
            int l=L+sizeTree, r=R+sizeTree;
            while(l<r){
                if(l&1){
                    segtree[l].push_back(item);
                    l++;
                }
                if(r&1){
                    --r;
                    segtree[r].push_back(item);
                }
                l>>=1; r>>=1;
            }
        }

        RollbackDSU dsu((int)edgeList.size());
        double best=1e100;
        bool found=false;

        function<void(int,int,int)> dfs = [&](int node,int lidx,int ridx){
            int checkpoint = (int)dsu.st.size();
            for(auto &e: segtree[node]){
                dsu.addEdge(e.u, e.v, e.p, e.q);
            }

            if(lidx>=N){
                dsu.rollback(checkpoint);
                return;
            }

            if(ridx-lidx==1){
                double repZ = ((double)H[lidx] + (double)H[lidx+1])*0.5;
                if(repZ > (double)overlap_low && repZ < (double)overlap_high){
                    int rw=dsu.find(westEdgeId);
                    int re=dsu.find(eastEdgeId);
                    if(rw==re){
                        double totP=dsu.sump[rw];
                        double totQ=dsu.sumq[rw];
                        double ZL = max((double)H[lidx], (double)overlap_low);
                        double ZR = min((double)H[lidx+1], (double)overlap_high);
                        double valL=totP*ZL + totQ;
                        double valR=totP*ZR + totQ;
                        double cand=min(valL,valR);
                        if(cand<0) cand=0;
                        if(cand<best){ best=cand; found=true; }
                    }
                }
            }else{
                int mid=(lidx+ridx)/2;
                dfs(node*2, lidx, mid);
                dfs(node*2+1, mid, ridx);
            }

            dsu.rollback(checkpoint);
        };

        dfs(1,0,sizeTree);

        if(!found){
            cout<<"impossible\n";
        } else {
            cout<<best<<"\n";
        }
    }
    return 0;
}