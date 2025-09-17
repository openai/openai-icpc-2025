#include <bits/stdc++.h>
using namespace std;

struct TriStep{
    int tri;
    int entry;
    double a;
    double b;
};

struct Triangle{
    int v[3];
    int neigh[3];
    int neighEdge[3];
    int lowLocal, midLocal, highLocal;
    int lowGlobal, midGlobal, highGlobal;
    double zLow, zMid, zHigh;
    double alphaLow, alphaHigh;
    bool active=false;
    bool modeLow=true;
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int tcase;
    if(!(cin>>tcase)) return 0;
    cout.setf(std::ios::fixed);
    cout<<setprecision(9);
    while(tcase--){
        int w,l;
        int n,m;
        cin>>w>>l>>n>>m;
        vector<double> x(n+1), y(n+1), z(n+1);
        int sw=-1,nw=-1,se=-1,ne=-1;
        for(int i=1;i<=n;i++){
            int xi,yi,zi;
            cin>>xi>>yi>>zi;
            x[i]=xi;y[i]=yi;z[i]=zi;
            if(xi==0 && yi==0) sw=i;
            if(xi==0 && yi==l) nw=i;
            if(xi==w && yi==0) se=i;
            if(xi==w && yi==l) ne=i;
        }

        vector<Triangle> tris(m);
        vector<vector<pair<int,int>>> incident(n+1);

        for(int t=0;t<m;t++){
            int a,b,c;
            cin>>a>>b>>c;
            tris[t].v[0]=a; tris[t].v[1]=b; tris[t].v[2]=c;
            for(int e=0;e<3;e++){
                tris[t].neigh[e]=-1;
                tris[t].neighEdge[e]=-1;
            }
            int ids[3]={a,b,c};
            sort(ids, ids+3, [&](int i1,int i2){return z[i1]<z[i2];});
            int low=ids[0], mid=ids[1], high=ids[2];
            tris[t].lowGlobal=low; tris[t].midGlobal=mid; tris[t].highGlobal=high;
            tris[t].zLow=z[low]; tris[t].zMid=z[mid]; tris[t].zHigh=z[high];
            int locLow=-1, locMid=-1, locHigh=-1;
            for(int k=0;k<3;k++){
                if(tris[t].v[k]==low) locLow=k;
                else if(tris[t].v[k]==mid) locMid=k;
                else if(tris[t].v[k]==high) locHigh=k;
            }
            tris[t].lowLocal=locLow; tris[t].midLocal=locMid; tris[t].highLocal=locHigh;
            incident[low].push_back({t,0});
            incident[mid].push_back({t,1});
            incident[high].push_back({t,2});

            double xl=x[low], yl=y[low];
            double xm=x[mid], ym=y[mid];
            double xh=x[high], yh=y[high];

            double dxHighLow=(xh - xl), dyHighLow=(yh - yl);
            double dxMidLow=(xm - xl), dyMidLow=(ym - yl);
            double Dlx = dxHighLow/(z[high]-z[low]) - dxMidLow/(z[mid]-z[low]);
            double Dly = dyHighLow/(z[high]-z[low]) - dyMidLow/(z[mid]-z[low]);
            double alphaLow = sqrt(Dlx*Dlx + Dly*Dly);
            tris[t].alphaLow=alphaLow;

            double dxmH = (xm - xh), dymH=(ym - yh);
            double dxlH = (xl - xh), dylH=(yl - yh);
            double Dhx = dxmH/(z[mid]-z[high]) - dxlH/(z[low]-z[high]);
            double Dhy = dymH/(z[mid]-z[high]) - dylH/(z[low]-z[high]);
            double alphaHigh = sqrt(Dhx*Dhx + Dhy*Dhy);
            tris[t].alphaHigh=alphaHigh;
        }

        unordered_map<long long, pair<int,int>> edgeMap;
        edgeMap.reserve(m*6+10);

        for(int tIdx=0;tIdx<m;tIdx++){
            for(int e=0;e<3;e++){
                int u=tris[tIdx].v[e];
                int vtx=tris[tIdx].v[(e+1)%3];
                int a=u,b=vtx;
                if(a>b) swap(a,b);
                long long key = ((long long)a<<17) | (long long)b;
                auto it=edgeMap.find(key);
                if(it==edgeMap.end()){
                    edgeMap[key]={tIdx,e};
                }else{
                    int ot=it->second.first;
                    int oe=it->second.second;
                    tris[tIdx].neigh[e]=ot;
                    tris[tIdx].neighEdge[e]=oe;
                    tris[ot].neigh[oe]=tIdx;
                    tris[ot].neighEdge[oe]=e;
                    edgeMap.erase(it);
                }
            }
        }

        int startTri=-1, startEntry=-1;
        for(int tIdx=0;tIdx<m;tIdx++){
            for(int e=0;e<3;e++){
                int u=tris[tIdx].v[e];
                int vtx=tris[tIdx].v[(e+1)%3];
                if((u==sw && vtx==nw) || (u==nw && vtx==sw)){
                    startTri=tIdx; startEntry=e; break;
                }
            }
            if(startTri!=-1) break;
        }

        vector<int> ord(n);
        iota(ord.begin(), ord.end(), 1);
        sort(ord.begin(), ord.end(), [&](int i, int j){return z[i] < z[j];});
        vector<double> Zv(n);
        for(int i=0;i<n;i++) Zv[i]=z[ord[i]];

        double westLow = min(z[sw], z[nw]);
        double westHigh = max(z[sw], z[nw]);
        double eastLow = min(z[se], z[ne]);
        double eastHigh = max(z[se], z[ne]);
        double Zlow = max(westLow, eastLow);
        double Zhigh = min(westHigh, eastHigh);

        vector<int> inPos(m,-1);
        vector<TriStep> chain;
        double a_tot=0.0, b_tot=0.0;
        bool endIsEast=false;
        bool found=false;
        double best=1e100;

        auto clearChain = [&](){
            while(!chain.empty()){
                TriStep st=chain.back();
                chain.pop_back();
                inPos[st.tri]=-1;
                a_tot-=st.a; b_tot-=st.b;
            }
            endIsEast=false;
        };

        auto getCoef = [&](int tIdx)->pair<double,double>{
            Triangle &tr=tris[tIdx];
            if(tr.modeLow){
                double a=tr.alphaLow;
                double b=-tr.alphaLow*tr.zLow;
                return {a,b};
            }else{
                double a=-tr.alphaHigh;
                double b=tr.alphaHigh*tr.zHigh;
                return {a,b};
            }
        };

        function<int(int,int)> getExit = [&](int tIdx, int entry)->int{
            Triangle &tr=tris[tIdx];
            int s = tr.modeLow? tr.lowLocal : tr.highLocal;
            int eA=s;
            int eB=(s+2)%3;
            if(entry==eA) return eB;
            else if(entry==eB) return eA;
            else return -1;
        };

        auto buildSuffix = [&](int startT, int startE){
            int currTri = startT;
            int currEntry = startE;
            while(currTri!=-1){
                if(inPos[currTri]!=-1){
                    endIsEast=false;
                    return;
                }
                Triangle &tr=tris[currTri];
                if(!tr.active) break;
                int s = tr.modeLow? tr.lowLocal : tr.highLocal;
                int eA=s;
                int eB=(s+2)%3;
                if(currEntry!=eA && currEntry!=eB){
                    break;
                }
                auto coef=getCoef(currTri);
                int idxNew= (int)chain.size();
                chain.push_back({currTri,currEntry,coef.first,coef.second});
                inPos[currTri]=idxNew;
                a_tot += coef.first;
                b_tot += coef.second;

                int exitEdge = (currEntry==eA)?eB:eA;
                int nextTri = tr.neigh[exitEdge];
                if(nextTri==-1){
                    int uu=tr.v[exitEdge];
                    int vv=tr.v[(exitEdge+1)%3];
                    if((uu==se && vv==ne) || (uu==ne && vv==se)) endIsEast=true;
                    else endIsEast=false;
                    currTri=-1;
                }else{
                    currEntry = tr.neighEdge[exitEdge];
                    currTri = nextTri;
                }
            }
        };

        for(int idxV=0; idxV<n; idxV++){
            int v=ord[idxV];
            int divergence=INT_MAX;
            for(auto &pr: incident[v]){
                int tIdx=pr.first;
                int role=pr.second;
                Triangle &tr=tris[tIdx];
                if(role==0){
                    tr.active=true;
                    tr.modeLow=true;
                }else if(role==1){
                    tr.modeLow=false;
                }else{
                    tr.active=false;
                }
                if(inPos[tIdx]!=-1){
                    divergence = min(divergence, inPos[tIdx]);
                }
            }

            if(idxV==n-1) break;
            double curZ=Zv[idxV];
            double nextZ=Zv[idxV+1];
            double rep = (curZ + nextZ)*0.5;
            bool westOk = (westLow < rep && rep < westHigh);
            bool eastOk = (eastLow < rep && rep < eastHigh);

            if(!westOk){
                if(!chain.empty()) clearChain();
            }else{
                if(chain.empty()){
                    if(tris[startTri].active) divergence=0;
                    else divergence=INT_MAX;
                }
                if(divergence!=INT_MAX){
                    while((int)chain.size()>divergence){
                        TriStep st=chain.back();
                        chain.pop_back();
                        inPos[st.tri]=-1;
                        a_tot-=st.a; b_tot-=st.b;
                    }
                    if(chain.empty()){
                        if(tris[startTri].active){
                            buildSuffix(startTri, startEntry);
                        }
                    }else{
                        TriStep &last=chain.back();
                        int exitE=getExit(last.tri, last.entry);
                        if(exitE==-1){
                            endIsEast=false;
                        }else{
                            int nextTri=tris[last.tri].neigh[exitE];
                            if(nextTri!=-1){
                                int entry=tris[last.tri].neighEdge[exitE];
                                buildSuffix(nextTri, entry);
                            }else{
                                int uu=tris[last.tri].v[exitE];
                                int vv=tris[last.tri].v[(exitE+1)%3];
                                if((uu==se && vv==ne) || (uu==ne && vv==se)) endIsEast=true;
                                else endIsEast=false;
                            }
                        }
                    }
                }
            }

            if(westOk && eastOk && !chain.empty() && endIsEast){
                if(curZ>=Zlow && curZ<=Zhigh){
                    double val=a_tot*curZ + b_tot;
                    if(val<best){best=val; found=true;}
                }
                if(nextZ>=Zlow && nextZ<=Zhigh){
                    double val=a_tot*nextZ + b_tot;
                    if(val<best){best=val; found=true;}
                }
            }
        }

        if(found){
            cout<<best<<"\n";
        }else{
            cout<<"impossible\n";
        }
    }
    return 0;
}