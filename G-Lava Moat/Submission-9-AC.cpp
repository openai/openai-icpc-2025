#include <bits/stdc++.h>
using namespace std;

struct Vert{double x,y,z;};
struct Edge{int u,v;};
struct Tri{int v[3]; int e[3]; int opp[3];};
struct ActEdge{int e1,e2; double a,b;};

struct DSURollback{
    vector<int> parent, sz;
    vector<double> A,B;
    struct Hist{int rv; int sz_ru; double A_ru; double B_ru;};
    vector<Hist> st;
    void init(int n){
        parent.resize(n);
        sz.resize(n);
        A.assign(n,0.0);
        B.assign(n,0.0);
        st.clear();
        for(int i=0;i<n;i++){
            parent[i]=i;
            sz[i]=1;
        }
    }
    int find(int x){
        while(parent[x]!=x) x=parent[x];
        return x;
    }
    bool unite(int u,int v,double a,double b){
        int ru=find(u), rv=find(v);
        if(ru==rv) return false;
        if(sz[ru]<sz[rv]) swap(ru,rv);
        st.push_back({rv,sz[ru],A[ru],B[ru]});
        parent[rv]=ru;
        sz[ru]+=sz[rv];
        A[ru] += A[rv] + a;
        B[ru] += B[rv] + b;
        return true;
    }
    void rollback(int target){
        while((int)st.size()>target){
            auto h=st.back(); st.pop_back();
            int rv=h.rv;
            int ru=parent[rv];
            parent[rv]=rv;
            sz[ru]=h.sz_ru;
            A[ru]=h.A_ru;
            B[ru]=h.B_ru;
        }
    }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t; 
    if(!(cin>>t)) return 0;
    cout.setf(std::ios::fixed); 
    cout<<setprecision(9);
    while(t--){
        int w,l; int n,m;
        cin>>w>>l>>n>>m;
        vector<Vert> vert(n+1);
        int swId=-1,seId=-1,nwId=-1,neId=-1;
        for(int i=1;i<=n;i++){
            int xi,yi,zi;
            cin>>xi>>yi>>zi;
            vert[i]={double(xi),double(yi),double(zi)};
            if(xi==0 && yi==0) swId=i;
            if(xi==w && yi==0) seId=i;
            if(xi==0 && yi==l) nwId=i;
            if(xi==w && yi==l) neId=i;
        }
        double minWest=min(vert[swId].z, vert[nwId].z);
        double maxWest=max(vert[swId].z, vert[nwId].z);
        double minEast=min(vert[seId].z, vert[neId].z);
        double maxEast=max(vert[seId].z, vert[neId].z);
        double overlapLow=max(minWest, minEast);
        double overlapHigh=min(maxWest, maxEast);

        unordered_map<long long,int> edgeMap;
        edgeMap.reserve(m*6);
        vector<Edge> edges;
        edges.reserve(m*3);
        vector<Tri> triangles;
        triangles.reserve(m);
        vector<vector<int>> incidentTris(n+1);

        auto getEdgeId = [&](int a,int b)->int{
            if(a>b) swap(a,b);
            long long key = (long long)a*60000LL + (long long)b;
            auto it=edgeMap.find(key);
            if(it!=edgeMap.end()) return it->second;
            int id=(int)edges.size();
            edges.push_back({a,b});
            edgeMap[key]=id;
            return id;
        };

        for(int ti=0;ti<m;ti++){
            int a,b,c1;
            cin>>a>>b>>c1;
            int e_ab=getEdgeId(a,b);
            int e_bc=getEdgeId(b,c1);
            int e_ca=getEdgeId(c1,a);
            Tri tr;
            tr.v[0]=a; tr.v[1]=b; tr.v[2]=c1;
            tr.e[0]=e_ab; tr.e[1]=e_bc; tr.e[2]=e_ca;
            tr.opp[0]=e_bc; tr.opp[1]=e_ca; tr.opp[2]=e_ab;
            int idx=(int)triangles.size();
            triangles.push_back(tr);
            incidentTris[a].push_back(idx);
            incidentTris[b].push_back(idx);
            incidentTris[c1].push_back(idx);
        }

        int E=(int)edges.size();
        int startEdgeId = getEdgeId(swId,nwId);
        int endEdgeId = getEdgeId(seId,neId);

        vector<int> order(n);
        iota(order.begin(), order.end(), 1);
        sort(order.begin(), order.end(), [&](int u,int v){ return vert[u].z < vert[v].z; });
        vector<int> pos(n+1);
        for(int i=0;i<n;i++) pos[order[i]]=i;

        int N=n;
        int sizeTree=1; while(sizeTree<N) sizeTree<<=1;
        vector<vector<ActEdge>> seg(2*sizeTree);

        auto addRange=[&](int L,int R,const ActEdge &ae){
            int l=L+sizeTree, r=R+sizeTree;
            while(l<r){
                if(l&1) seg[l++].push_back(ae);
                if(r&1) seg[--r].push_back(ae);
                l>>=1; r>>=1;
            }
        };

        for(const auto &tr: triangles){
            int ids[3]={tr.v[0],tr.v[1],tr.v[2]};
            int ordIdx[3]={0,1,2};
            sort(ordIdx, ordIdx+3, [&](int ii,int jj){ return vert[ ids[ii] ].z < vert[ ids[jj] ].z; });
            int lo=ids[ordIdx[0]];
            int mi=ids[ordIdx[1]];
            int hi=ids[ordIdx[2]];
            int plo=pos[lo], pmi=pos[mi], phi=pos[hi];

            int L=plo+1, R=pmi;
            if(L<R){
                int sing=lo; int q=mi; int r=hi;
                const Vert &S=vert[sing], &Q=vert[q], &RR=vert[r];
                double dx1=Q.x-S.x, dy1=Q.y-S.y;
                double dx2=RR.x-S.x, dy2=RR.y-S.y;
                double denom1=Q.z-S.z, denom2=RR.z-S.z;
                double kx=dx1/denom1 - dx2/denom2;
                double ky=dy1/denom1 - dy2/denom2;
                double normK=sqrt(kx*kx+ky*ky);
                double aCoeff=+normK;
                double bCoeff=-S.z*normK;
                int e1=-1,e2=-1;
                for(int j=0;j<3;j++){
                    int eid=tr.e[j];
                    auto &ed=edges[eid];
                    if(ed.u==sing||ed.v==sing){
                        if(e1==-1) e1=eid; else e2=eid;
                    }
                }
                ActEdge ae{e1,e2,aCoeff,bCoeff};
                addRange(L,R,ae);
            }

            L=pmi+1; R=phi;
            if(L<R){
                int sing=hi; int q=lo; int r=mi;
                const Vert &S=vert[sing], &Q=vert[q], &RR=vert[r];
                double dx1=Q.x-S.x, dy1=Q.y-S.y;
                double dx2=RR.x-S.x, dy2=RR.y-S.y;
                double denom1=Q.z-S.z, denom2=RR.z-S.z;
                double kx=dx1/denom1 - dx2/denom2;
                double ky=dy1/denom1 - dy2/denom2;
                double normK=sqrt(kx*kx+ky*ky);
                double aCoeff=-normK;
                double bCoeff=+S.z*normK;
                int e1=-1,e2=-1;
                for(int j=0;j<3;j++){
                    int eid=tr.e[j];
                    auto &ed=edges[eid];
                    if(ed.u==sing||ed.v==sing){
                        if(e1==-1) e1=eid; else e2=eid;
                    }
                }
                ActEdge ae{e1,e2,aCoeff,bCoeff};
                addRange(L,R,ae);
            }
        }

        DSURollback dsu;
        dsu.init(E);
        double ans=numeric_limits<double>::infinity();

        function<void(int,int,int)> dfs = [&](int idx,int l,int r){
            int before=(int)dsu.st.size();
            for(auto &ae: seg[idx]){
                dsu.unite(ae.e1, ae.e2, ae.a, ae.b);
            }

            if(r-l==1){
                int leaf=l;
                if(leaf < N){
                    int vid=order[leaf];
                    double c = vert[vid].z;
                    if(!(c < overlapLow || c > overlapHigh)){
                        int rootStart=dsu.find(startEdgeId);
                        int rootEnd=dsu.find(endEdgeId);
                        double lenStart=dsu.A[rootStart]*c + dsu.B[rootStart];
                        double lenEnd=dsu.A[rootEnd]*c + dsu.B[rootEnd];

                        double best=numeric_limits<double>::infinity();

                        if(rootStart==rootEnd){
                            best=lenStart;
                        }

                        double minStartCost=numeric_limits<double>::infinity();
                        double minEndCost=numeric_limits<double>::infinity();
                        if(vid==swId || vid==nwId) minStartCost=0.0;
                        if(vid==seId || vid==neId) minEndCost=0.0;

                        for(int tIdx: incidentTris[vid]){
                            const Tri &tr=triangles[tIdx];
                            int idxv;
                            if(tr.v[0]==vid) idxv=0;
                            else if(tr.v[1]==vid) idxv=1;
                            else idxv=2;
                            int a=tr.v[(idxv+1)%3];
                            int b=tr.v[(idxv+2)%3];
                            double za=vert[a].z, zb=vert[b].z;
                            double da=za - c, db=zb - c;
                            if(!((da>0 && db<0) || (da<0 && db>0))) continue;
                            int eOpp=tr.opp[idxv];
                            int rootOpp=dsu.find(eOpp);
                            double lenOppRoot=dsu.A[rootOpp]*c + dsu.B[rootOpp];
                            const Vert &Av=vert[a]; const Vert &Bv=vert[b]; const Vert &Vv=vert[vid];
                            double tpar=(c - Av.z) / (Bv.z - Av.z);
                            double px=Av.x + tpar*(Bv.x - Av.x);
                            double py=Av.y + tpar*(Bv.y - Av.y);
                            double segLen=hypot(px - Vv.x, py - Vv.y);

                            if(rootOpp==rootStart){
                                double cost = ((eOpp==startEdgeId)?0.0:lenStart) + segLen;
                                if(cost<minStartCost) minStartCost=cost;
                            }
                            if(rootOpp==rootEnd){
                                double cost = ((eOpp==endEdgeId)?0.0:lenEnd) + segLen;
                                if(cost<minEndCost) minEndCost=cost;
                            }
                        }

                        if(minStartCost<numeric_limits<double>::infinity() && minEndCost<numeric_limits<double>::infinity()){
                            double tot=minStartCost+minEndCost;
                            if(tot<best) best=tot;
                        }

                        if(best<ans) ans=best;
                    }
                }
            }else{
                int mid=(l+r)/2;
                dfs(idx*2,l,mid);
                dfs(idx*2+1,mid,r);
            }

            dsu.rollback(before);
        };

        dfs(1,0,sizeTree);

        if(std::isinf(ans)){
            cout<<"impossible\n";
        }else{
            cout<<ans<<"\n";
        }
    }

    return 0;
}