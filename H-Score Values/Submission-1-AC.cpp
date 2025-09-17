#include <bits/stdc++.h>
using namespace std;

using ll = long long;

ll gcdll(ll a, ll b){ return b?gcdll(b,a%b):a; }

vector<int> countDigitsLL(long long x){
    vector<int> cnt(10,0);
    if(x==0){ cnt[0]=1; return cnt; }
    while(x>0){ cnt[x%10]++; x/=10; }
    return cnt;
}

vector<int> toDigitsPadded(long long x, int len){
    string s = to_string(x);
    if((int)s.size() < len) s = string(len - s.size(), '0') + s;
    vector<int> v(len);
    for(int i=0;i<len;i++) v[i] = s[i]-'0';
    return v;
}

// Global for DP
int Gmod; // modulus
vector<int> digM, digL; int LEN;
long long w_digits[10];
static long long memo[20][1001][2][2][2];
static unsigned char vis[20][1001][2][2][2];
const long long INFLL = (1LL<<60);

long long dfsDP(int pos, int rem, int tightL, int tightU, int started){
    if(pos==LEN){
        if(rem % Gmod == 0){
            if(!started) return w_digits[0];
            else return 0LL;
        }
        return -INFLL/4;
    }
    long long &res = memo[pos][rem][tightL][tightU][started];
    if(vis[pos][rem][tightL][tightU][started]) return res;
    vis[pos][rem][tightL][tightU][started] = 1;
    long long best = -INFLL/4;
    int lowdigit = digL[pos];
    int updigit = digM[pos];
    int lo = tightL? lowdigit : 0;
    int hi = tightU? updigit : 9;
    if(lo>hi){ res = -INFLL/4; return res; }
    for(int d=lo; d<=hi; ++d){
        int newStarted = started || (d!=0);
        long long add = newStarted ? w_digits[d] : 0LL;
        int newRem = ((long long)rem*10 + d) % Gmod;
        int newTightL = tightL && (d==lowdigit);
        int newTightU = tightU && (d==updigit);
        long long val = dfsDP(pos+1, newRem, newTightL, newTightU, newStarted);
        if(val > -INFLL/8){
            long long cand = add + val;
            if(cand > best) best = cand;
        }
    }
    res = best;
    return res;
}

long long dpMaxWeight(long long L, long long M, int G, const vector<int>& digM_in, const vector<int>& digL_in, int len, const vector<int>& w){
    if(L>M) return -INFLL/4;
    Gmod = max(1, G);
    LEN = len; digM = digM_in; digL = digL_in;
    for(int i=0;i<10;i++) w_digits[i] = w[i];
    memset(vis,0,sizeof(vis));
    long long ans = dfsDP(0,0,1,1,0);
    return ans;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    long long m; int n;
    if(!(cin>>m>>n)) return 0;
    vector<int> p(n);
    for(int i=0;i<n;i++) cin>>p[i];

    ll g=0; for(int x: p) g = gcdll(g,(ll)x);
    vector<int> q(n); for(int i=0;i<n;i++) q[i] = p[i]/g;

    int a_min = *min_element(q.begin(), q.end());
    int N = a_min;
    vector<long long> dist(N, (long long)4e18);
    using P=pair<long long,int>;
    priority_queue<P, vector<P>, greater<P>> pq;
    dist[0]=0; pq.push({0,0});
    while(!pq.empty()){
        auto [d,u]=pq.top(); pq.pop();
        if(d!=dist[u]) continue;
        for(int w: q){
            int v = (u+w)%N;
            long long nd = d + w;
            if(nd < dist[v]){
                dist[v] = nd;
                pq.push({nd,v});
            }
        }
    }
    long long Bprime = 0;
    for(long long dv: dist) if(dv>Bprime) Bprime = dv;

    long long Ymax = m / g;

    vector<long long> ans(9, 0);

    long long limitSmall = min(Ymax, Bprime-1);
    if(Bprime==0) limitSmall = -1;
    for(long long y=0; y<=limitSmall; ++y){
        if(y >= dist[(int)(y%N)]){
            long long x = y * g;
            vector<int> cnt = countDigitsLL(x);
            for(int d=0; d<=8; ++d){
                if(d==6){ long long c = (long long)cnt[6] + cnt[9]; if(c>ans[d]) ans[d]=c; }
                else { long long c = cnt[d]; if(c>ans[d]) ans[d]=c; }
            }
        }
    }

    long long Lx = g * Bprime;

    string sM = to_string(m); int len = sM.size(); 
    vector<int> digM_all(len); for(int i=0;i<len;i++) digM_all[i]=sM[i]-'0';
    vector<int> digL_all = toDigitsPadded(Lx, len);

    if(Ymax >= Bprime){
        for(int d=0; d<=8; ++d){
            vector<int> w(10,0);
            if(d==6){ w[6]=1; w[9]=1; }
            else { w[d]=1; }
            long long best = dpMaxWeight(Lx, m, (int)g, digM_all, digL_all, len, w);
            if(best > ans[d]) ans[d] = best;
        }
    }

    // Add m itself due to cap
    vector<int> cntm = countDigitsLL(m);
    for(int d=0; d<=8; ++d){
        if(d==6){ long long c = (long long)cntm[6] + cntm[9]; if(c>ans[d]) ans[d]=c; }
        else { long long c = cntm[d]; if(c>ans[d]) ans[d]=c; }
    }

    for(int d=0; d<=8; ++d){
        if(ans[d] > 0) cout<<d<<" "<<ans[d]<<"\n";
    }

    return 0;
}