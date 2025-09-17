#include <bits/stdc++.h>
using namespace std;

struct Change{
    int rv;
    int ru;
    int oldSizeRu;
    double oldKru;
    double oldBru;
    bool merged;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    if(!(cin >> t)) return 0;
    cout.setf(std::ios::fixed);
    cout << setprecision(9);
    while(t--){
        int w,l;
        int n,m;
        cin >> w >> l >> n >> m;
        vector<double> x(n+1), y(n+1), z(n+1);
        int sw=-1,se=-1,nw=-1,ne=-1;
        for(int i=1;i<=n;i++){
            int xi,yi,zi;
            cin >> xi >> yi >> zi;
            x[i]=xi; y[i]=yi; z[i]=zi;
            if(xi==0 && yi==0) sw=i;
            if(xi==w && yi==0) se=i;
            if(xi==w && yi==l) ne=i;
            if(xi==0 && yi==l) nw=i;
        }

        unordered_map<long long,int> edgeId;
        edgeId.reserve(3*m*2);
        vector<pair<int,int>> edgeEndpoints;
        edgeEndpoints.reserve(3*m);

        auto addEdge = [&](int a,int b)->int{
            if(a>b) swap(a,b);
            long long key = (long long)a*((long long)n+1LL) + (long long)b;
            auto it=edgeId.find(key);
            if(it==edgeId.end()){
                int id=(int)edgeEndpoints.size();
                edgeEndpoints.push_back({a,b});
                edgeId[key]=id;
                return id;
            }
            return it->second;
        };
        auto lookupEdge = [&](int a,int b)->int{
            if(a>b) swap(a,b);
            long long key = (long long)a*((long long)n+1LL) + (long long)b;
            auto it=edgeId.find(key);
            if(it==edgeId.end()) return -1;
            return it->second;
        };

        struct Tri{int v[3];};
        vector<Tri> triangles;
        triangles.reserve(m);

        for(int i=0;i<m;i++){
            int a,b,c;
            cin >> a >> b >> c;
            Tri tr; tr.v[0]=a; tr.v[1]=b; tr.v[2]=c;
            triangles.push_back(tr);
            addEdge(a,b);
            addEdge(b,c);
            addEdge(c,a);
        }

        int N=n;
        vector<int> idxSorted(N);
        iota(idxSorted.begin(), idxSorted.end(), 1);
        sort(idxSorted.begin(), idxSorted.end(), [&](int u,int v){return z[u] < z[v];});
        vector<int> ord(n+1);
        vector<double> hSorted(N);
        for(int i=0;i<N;i++){
            ord[idxSorted[i]]=i;
            hSorted[i]=z[idxSorted[i]];
        }

        int westEdgeId = lookupEdge(sw,nw);
        int eastEdgeId = lookupEdge(se,ne);
        if(westEdgeId<0 || eastEdgeId<0){
            cout << "impossible\n";
            continue;
        }

        vector<int> gu, gv;
        vector<double> gK, gB;
        vector<int> gStart, gEnd;
        gu.reserve(2*m);
        gv.reserve(2*m);
        gK.reserve(2*m);
        gB.reserve(2*m);
        gStart.reserve(2*m);
        gEnd.reserve(2*m);

        for(auto &tr: triangles){
            int a=tr.v[0], b=tr.v[1], c=tr.v[2];
            int v0=a, v1=b, v2=c;
            if(z[v1] < z[v0]) swap(v0,v1);
            if(z[v2] < z[v1]) swap(v1,v2);
            if(z[v1] < z[v0]) swap(v0,v1);
            int lo=v0, mi=v1, hi=v2;

            double zlo=z[lo], zmi=z[mi], zhi=z[hi];
            double xlo=x[lo], ylo=y[lo];
            double xmi=x[mi], ymi=y[mi];
            double xhi=x[hi], yhi=y[hi];

            int plo=ord[lo], pmi=ord[mi], phi=ord[hi];

            int start1=plo+1, end1=pmi+1;
            int e_lomi=lookupEdge(lo,mi);
            int e_lohi=lookupEdge(lo,hi);
            if(e_lomi>=0 && e_lohi>=0 && start1<end1){
                double dx = (xmi - xlo)/(zmi - zlo) - (xhi - xlo)/(zhi - zlo);
                double dy = (ymi - ylo)/(zmi - zlo) - (yhi - ylo)/(zhi - zlo);
                double K1 = sqrt(dx*dx + dy*dy);
                double B1 = -K1*zlo;
                gu.push_back(e_lomi);
                gv.push_back(e_lohi);
                gK.push_back(K1);
                gB.push_back(B1);
                gStart.push_back(start1);
                gEnd.push_back(end1);
            }

            int start2=pmi+1, end2=phi+1;
            int e_mihi=lookupEdge(mi,hi);
            if(e_mihi>=0 && e_lohi>=0 && start2<end2){
                double dx = (xlo - xhi)/(zhi - zlo) - (xmi - xhi)/(zhi - zmi);
                double dy = (ylo - yhi)/(zhi - zlo) - (ymi - yhi)/(zhi - zmi);
                double K2 = sqrt(dx*dx + dy*dy);
                double Kcoef = -K2;
                double B2 = K2*zhi;
                gu.push_back(e_lohi);
                gv.push_back(e_mihi);
                gK.push_back(Kcoef);
                gB.push_back(B2);
                gStart.push_back(start2);
                gEnd.push_back(end2);
            }
        }

        int NStates = N + 1;
        int Msize=1;
        while(Msize < NStates) Msize<<=1;
        vector<vector<int>> seg(2*Msize);

        int gNum=gu.size();
        for(int idx=0; idx<gNum; idx++){
            int l=gStart[idx];
            int r=gEnd[idx];
            if(l<0) l=0;
            if(r> NStates) r=NStates;
            if(l>=r) continue;
            int L=max(l,0), R=min(r,NStates);
            if(L>=R) continue;
            int ll=L+Msize, rr=R+Msize;
            while(ll<rr){
                if(ll&1) seg[ll++].push_back(idx);
                if(rr&1) seg[--rr].push_back(idx);
                ll>>=1; rr>>=1;
            }
        }

        int V=edgeEndpoints.size();
        vector<int> parent(V), sz(V);
        vector<double> kSum(V,0.0), bSum(V,0.0);
        for(int i=0;i<V;i++){parent[i]=i; sz[i]=1;}
        vector<Change> st;
        st.reserve(gNum);

        function<int(int)> findRoot = [&](int x){
            while(parent[x]!=x) x=parent[x];
            return x;
        };

        auto addG = [&](int idx){
            int u=gu[idx], v=gv[idx];
            double K=gK[idx], B=gB[idx];
            int ru=findRoot(u), rv=findRoot(v);
            if(ru==rv){
                Change ch{ru, ru, 0, kSum[ru], bSum[ru], false};
                st.push_back(ch);
                kSum[ru] += K;
                bSum[ru] += B;
            }else{
                if(sz[ru] < sz[rv]) swap(ru,rv);
                Change ch{rv, ru, sz[ru], kSum[ru], bSum[ru], true};
                st.push_back(ch);
                parent[rv]=ru;
                sz[ru] += sz[rv];
                kSum[ru] += kSum[rv] + K;
                bSum[ru] += bSum[rv] + B;
            }
        };

        auto rollback = [&](int checkpoint){
            while((int)st.size() > checkpoint){
                Change ch=st.back();
                st.pop_back();
                if(!ch.merged){
                    int r=ch.rv;
                    kSum[r]=ch.oldKru;
                    bSum[r]=ch.oldBru;
                }else{
                    int rv=ch.rv, ru=ch.ru;
                    parent[rv]=rv;
                    sz[ru]=ch.oldSizeRu;
                    kSum[ru]=ch.oldKru;
                    bSum[ru]=ch.oldBru;
                }
            }
        };

        vector<double> stateK(NStates,0.0), stateB(NStates,0.0);
        vector<char> stateConnected(NStates,0);

        function<void(int,int,int)> dfs = [&](int node,int l,int r){
            int checkpoint = (int)st.size();
            for(int idx: seg[node]) addG(idx);
            if(r-l==1){
                int k=l;
                if(k < NStates){
                    int rw=findRoot(westEdgeId);
                    int re=findRoot(eastEdgeId);
                    stateK[k]=kSum[rw];
                    stateB[k]=bSum[rw];
                    stateConnected[k]=(rw==re);
                }
            }else{
                int mid=(l+r)/2;
                dfs(node*2,l,mid);
                dfs(node*2+1,mid,r);
            }
            rollback(checkpoint);
        };

        dfs(1,0,Msize);

        auto isActiveEdge = [&](int a,int b,int k)->bool{
            int oa=ord[a], ob=ord[b];
            if(oa>ob) swap(oa,ob);
            return (oa < k && k <= ob);
        };

        double ans=1e100;
        for(int p=0;p<N;p++){
            double cvalue=hSorted[p];
            int kLeft=p;
            if(kLeft < NStates){
                if(isActiveEdge(sw,nw,kLeft) && isActiveEdge(se,ne,kLeft) && stateConnected[kLeft]){
                    double dist=stateK[kLeft]*cvalue + stateB[kLeft];
                    if(dist<ans) ans=dist;
                }
            }
            int kRight=p+1;
            if(kRight < NStates){
                if(isActiveEdge(sw,nw,kRight) && isActiveEdge(se,ne,kRight) && stateConnected[kRight]){
                    double dist=stateK[kRight]*cvalue + stateB[kRight];
                    if(dist<ans) ans=dist;
                }
            }
        }

        if(ans==1e100){
            cout << "impossible\n";
        }else{
            cout << setprecision(9) << ans << "\n";
        }
    }
    return 0;
}