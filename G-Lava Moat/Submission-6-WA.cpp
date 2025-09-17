#include <bits/stdc++.h>
using namespace std;

struct Vertex{double x,y,z;};
struct TriData{int v[3]; double a1,b1,a2,b2;};
struct Edge{int u,v; int tri1; int tri2;};

struct Node{
    int id;
    int prio;
    Node* left;
    Node* right;
    Node* parent;
    bool rev;
    double valA,valB;
    double sumA,sumB;
    int sz;
    bool isCycle;
};

static inline int getSize(Node* t){return t?t->sz:0;}
static inline void upd(Node* t){
    if(!t) return;
    t->sz=1;
    t->sumA=t->valA;
    t->sumB=t->valB;
    if(t->left){
        t->sz += t->left->sz;
        t->sumA += t->left->sumA;
        t->sumB += t->left->sumB;
    }
    if(t->right){
        t->sz += t->right->sz;
        t->sumA += t->right->sumA;
        t->sumB += t->right->sumB;
    }
}
static inline void push(Node* t){
    if(t && t->rev){
        t->rev=false;
        Node* tmp=t->left; t->left=t->right; t->right=tmp;
        if(t->left) t->left->rev^=1;
        if(t->right) t->right->rev^=1;
    }
}

static Node* mergeTreap(Node* a, Node* b){
    if(!a){ if(b) b->parent=nullptr; return b;}
    if(!b){ a->parent=nullptr; return a;}
    if(a->prio > b->prio){
        push(a);
        Node* mergedRight = mergeTreap(a->right,b);
        a->right = mergedRight;
        if(a->right) a->right->parent=a;
        upd(a);
        return a;
    }else{
        push(b);
        Node* mergedLeft = mergeTreap(a,b->left);
        b->left = mergedLeft;
        if(b->left) b->left->parent=b;
        upd(b);
        return b;
    }
}

static pair<Node*,Node*> splitTreap(Node* t, int k){
    if(!t) return {nullptr,nullptr};
    push(t);
    int lsize=getSize(t->left);
    if(k<=lsize){
        auto res=splitTreap(t->left,k);
        t->left=res.second;
        if(t->left) t->left->parent=t;
        upd(t);
        t->parent=nullptr;
        if(res.first) res.first->parent=nullptr;
        return {res.first,t};
    }else{
        auto res=splitTreap(t->right,k-lsize-1);
        t->right=res.first;
        if(t->right) t->right->parent=t;
        upd(t);
        if(res.second) res.second->parent=nullptr;
        t->parent=nullptr;
        return {t,res.second};
    }
}

static Node* getRoot(Node* x){
    while(x->parent) x=x->parent;
    return x;
}

static int getIndex(Node* x){
    vector<Node*> path;
    Node* cur=x;
    while(cur){
        path.push_back(cur);
        cur=cur->parent;
    }
    for(int i=(int)path.size()-1;i>=0;--i) push(path[i]);
    int idx=getSize(x->left)+1;
    cur=x;
    while(cur->parent){
        Node* p=cur->parent;
        if(cur==p->right){
            idx += getSize(p->left)+1;
        }
        cur=p;
    }
    return idx;
}

struct Solver{
    int w,l;
    int n,m;
    vector<Vertex> verts;
    vector<TriData> tris;
    vector<Edge> edges;
    unordered_map<long long,int> edgeMap;
    vector<vector<int>> adjEdges;
    vector<vector<int>> trisOfVertex;
    vector<char> inserted;
    int sw,nw,se,ne;
    int westEdgeId,eastEdgeId;
    int startTri,endTri;
    int N;
    int LEFT_ID, RIGHT_ID;
    vector<Node*> nodes;
    vector<int> deg;
    vector<array<int,2>> neigh;
    mt19937 rng;

    Solver():rng(712367123){}

    long long makeKey(int a,int b){
        if(a>b) swap(a,b);
        return ((long long)a<<17) | (long long)b;
    }

    void addEdgeToMap(int u,int v,int triId){
        long long key=makeKey(u,v);
        auto it=edgeMap.find(key);
        if(it==edgeMap.end()){
            int id=edges.size();
            edges.push_back({u,v,triId,-1});
            edgeMap[key]=id;
            adjEdges[u].push_back(id);
            adjEdges[v].push_back(id);
        }else{
            int id=it->second;
            edges[id].tri2=triId;
        }
    }

    void initTest(){
        verts.assign(n,{});
        tris.assign(m,{});
        edges.clear();
        edgeMap.clear();
        adjEdges.assign(n,{});
        trisOfVertex.assign(n,{});
        inserted.assign(n,0);
    }

    void buildTreapNodes(){
        N=m+2; LEFT_ID=m; RIGHT_ID=m+1;
        nodes.assign(N,nullptr);
        deg.assign(N,0);
        neigh.assign(N,{});
        vector<int> prios(N);
        iota(prios.begin(), prios.end(), 0);
        shuffle(prios.begin(), prios.end(), rng);
        for(int i=0;i<N;i++){
            Node* nd=new Node();
            nd->id=i;
            nd->prio=prios[i];
            nd->left=nd->right=nd->parent=nullptr;
            nd->rev=false;
            nd->valA=nd->valB=0;
            nd->sumA=nd->sumB=0;
            nd->sz=1;
            nd->isCycle=false;
            nodes[i]=nd;
        }
    }

    void addNeighbor(int u,int v){
        int d=deg[u];
        neigh[u][d]=v;
        deg[u]=d+1;
    }
    void removeNeighbor(int u,int v){
        int d=deg[u];
        if(d==0) return;
        if(neigh[u][0]==v){
            if(d==2) neigh[u][0]=neigh[u][1];
            deg[u]--;
        }else if(d==2 && neigh[u][1]==v){
            deg[u]--;
        }
    }

    void linkTreap(int u,int v){
        Node* nu=nodes[u];
        Node* nv=nodes[v];
        Node* ru=getRoot(nu);
        Node* rv=getRoot(nv);
        if(ru!=rv){
            int su=ru->sz;
            int sv=rv->sz;
            int idxu=getIndex(nu);
            int idxv=getIndex(nv);
            if(idxu!=su){
                ru->rev^=1;
            }
            if(idxv!=1){
                rv->rev^=1;
            }
            Node* merged=mergeTreap(ru,rv);
            if(merged) merged->parent=nullptr;
            merged->isCycle=false;
        }else{
            ru->isCycle=true;
        }
    }

    void cutTreap(int u,int v){
        Node* nu=nodes[u];
        Node* nv=nodes[v];
        Node* root=getRoot(nu);
        int idxu=getIndex(nu);
        int idxv=getIndex(nv);
        int nsize=root->sz;
        if(!root->isCycle){
            int pos=min(idxu,idxv);
            auto pr=splitTreap(root,pos);
            Node* left=pr.first;
            Node* right=pr.second;
            if(left){left->isCycle=false; left->parent=nullptr;}
            if(right){right->isCycle=false; right->parent=nullptr;}
        }else{
            if((idxu==1 && idxv==nsize) || (idxv==1 && idxu==nsize)){
                root->isCycle=false;
            }else{
                int pos=min(idxu,idxv);
                auto pr=splitTreap(root,pos);
                Node* left=pr.first;
                Node* right=pr.second;
                Node* merged=mergeTreap(right,left);
                if(merged){merged->parent=nullptr; merged->isCycle=false;}
            }
        }
    }

    void addGraphEdge(int u,int v){
        addNeighbor(u,v);
        addNeighbor(v,u);
        linkTreap(u,v);
    }
    void removeGraphEdge(int u,int v){
        removeNeighbor(u,v);
        removeNeighbor(v,u);
        cutTreap(u,v);
    }

    void handleAddCut(int eId){
        Edge &e=edges[eId];
        if(eId==westEdgeId){
            int t=e.tri1!=-1?e.tri1:e.tri2;
            addGraphEdge(LEFT_ID,t);
        }else if(eId==eastEdgeId){
            int t=e.tri1!=-1?e.tri1:e.tri2;
            addGraphEdge(RIGHT_ID,t);
        }else{
            if(e.tri2!=-1){
                addGraphEdge(e.tri1,e.tri2);
            }
        }
    }

    void handleRemoveCut(int eId){
        Edge &e=edges[eId];
        if(eId==westEdgeId){
            int t=e.tri1!=-1?e.tri1:e.tri2;
            removeGraphEdge(LEFT_ID,t);
        }else if(eId==eastEdgeId){
            int t=e.tri1!=-1?e.tri1:e.tri2;
            removeGraphEdge(RIGHT_ID,t);
        }else{
            if(e.tri2!=-1){
                removeGraphEdge(e.tri1,e.tri2);
            }
        }
    }

    void updateWeight(int triId,double newA,double newB){
        Node* nd=nodes[triId];
        nd->valA=newA;
        nd->valB=newB;
        Node* cur=nd;
        while(cur){
            upd(cur);
            cur=cur->parent;
        }
    }

    double solveOne(){
        long long keyW=makeKey(sw,nw);
        westEdgeId=edgeMap[keyW];
        long long keyE=makeKey(se,ne);
        eastEdgeId=edgeMap[keyE];
        startTri = edges[westEdgeId].tri1!=-1?edges[westEdgeId].tri1:edges[westEdgeId].tri2;
        endTri = edges[eastEdgeId].tri1!=-1?edges[eastEdgeId].tri1:edges[eastEdgeId].tri2;

        buildTreapNodes();

        vector<int> order(n);
        iota(order.begin(),order.end(),0);
        sort(order.begin(),order.end(),[&](int a,int b){return verts[a].z<verts[b].z;});

        double ans=numeric_limits<double>::infinity();

        for(int idx=0; idx<n; ++idx){
            int v=order[idx];
            double h=verts[v].z;

            // evaluate before insertion
            {
                Node* rootL=getRoot(nodes[LEFT_ID]);
                Node* rootR=getRoot(nodes[RIGHT_ID]);
                if(rootL==rootR){
                    double cand = rootL->sumA * h + rootL->sumB;
                    if(cand<ans) ans=cand;
                }
            }

            vector<int> rem, add;
            for(int eId: adjEdges[v]){
                Edge &e=edges[eId];
                int u = (e.u==v)?e.v:e.u;
                if(inserted[u]) rem.push_back(eId);
                else add.push_back(eId);
            }

            for(int eId: rem) handleRemoveCut(eId);
            for(int eId: add) handleAddCut(eId);

            inserted[v]=1;
            for(int tId: trisOfVertex[v]){
                TriData &tr=tris[tId];
                int c = (int)inserted[tr.v[0]] + (int)inserted[tr.v[1]] + (int)inserted[tr.v[2]];
                double na=0, nb=0;
                if(c==1){na=tr.a1; nb=tr.b1;}
                else if(c==2){na=tr.a2; nb=tr.b2;}
                updateWeight(tId,na,nb);
            }

            // evaluate after insertion
            {
                Node* rootL=getRoot(nodes[LEFT_ID]);
                Node* rootR=getRoot(nodes[RIGHT_ID]);
                if(rootL==rootR){
                    double cand = rootL->sumA * h + rootL->sumB;
                    if(cand<ans) ans=cand;
                }
            }
        }

        return ans;
    }

    void run(){
        ios::sync_with_stdio(false);
        cin.tie(nullptr);
        int tcase;
        if(!(cin>>tcase)) return;
        cout.setf(std::ios::fixed);
        cout<<setprecision(9);
        while(tcase--){
            cin>>w>>l>>n>>m;
            initTest();
            sw=nw=se=ne=-1;
            for(int i=0;i<n;i++){
                long long x,y,z;
                cin>>x>>y>>z;
                verts[i]={double(x),double(y),double(z)};
                if(x==0 && y==0) sw=i;
                if(x==0 && y==l) nw=i;
                if(x==w && y==0) se=i;
                if(x==w && y==l) ne=i;
            }
            for(int ti=0;ti<m;ti++){
                int a,b,c;
                cin>>a>>b>>c;
                a--;b--;c--;
                tris[ti].v[0]=a; tris[ti].v[1]=b; tris[ti].v[2]=c;
                trisOfVertex[a].push_back(ti);
                trisOfVertex[b].push_back(ti);
                trisOfVertex[c].push_back(ti);
                addEdgeToMap(a,b,ti);
                addEdgeToMap(b,c,ti);
                addEdgeToMap(c,a,ti);
            }

            for(int ti=0;ti<m;ti++){
                TriData &tr=tris[ti];
                int ids[3] = {tr.v[0],tr.v[1],tr.v[2]};
                int ord[3] = {0,1,2};
                sort(ord,ord+3,[&](int ii,int jj){return verts[ids[ii]].z < verts[ids[jj]].z;});
                int v0=ids[ord[0]], v1=ids[ord[1]], v2=ids[ord[2]];
                double z0=verts[v0].z, z1=verts[v1].z, z2=verts[v2].z;
                double x0=verts[v0].x, y0=verts[v0].y;
                double x1=verts[v1].x, y1=verts[v1].y;
                double x2=verts[v2].x, y2=verts[v2].y;
                double tparam=(z1-z0)/(z2-z0);
                double px=x0 + tparam*(x2-x0);
                double py=y0 + tparam*(y2-y0);
                double dist=hypot(px-x1,py-y1);
                double a1=dist/(z1-z0);
                double b1=-a1*z0;
                double a2=-dist/(z2-z1);
                double b2=dist*z2/(z2-z1);
                tr.a1=a1; tr.b1=b1; tr.a2=a2; tr.b2=b2;
            }

            double ans=solveOne();
            if(!isfinite(ans)){
                cout<<"impossible"<<"\n";
            }else{
                cout<<ans<<"\n";
            }
            for(Node* nd: nodes) delete nd;
        }
    }
};

int main(){
    Solver solver;
    solver.run();
    return 0;
}