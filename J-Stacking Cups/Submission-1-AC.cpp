#include <bits/stdc++.h>
using namespace std;
using ll = long long;

static inline long long isqrtll(long long x){
    long long r = sqrt((long double)x);
    while((r+1)*(r+1) <= x) r++;
    while(r*r > x) r--;
    return r;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    long long n, h;
    if(!(cin >> n >> h)) return 0;

    long long minH = 2*n - 1;
    long long maxH = n*n;

    if(h < minH || h > maxH || (n>=3 && h==maxH-2)){
        cout << "impossible\n";
        return 0;
    }

    if(h == minH){
        for(long long i=n;i>=1;i--){
            if(i!=n) cout << ' ';
            cout << (2*i-1);
        }
        cout << '\n';
        return 0;
    }

    vector<int> order; order.reserve(n);

    for(int b=1;b<=n;b++){
        long long t = h - (n - b);
        if(t < 0 || t > 1LL*b*b) continue;
        long long a = isqrtll(t);
        if((a & 1) != (t & 1)) a--;
        if(a < 0) continue;
        if(t <= a*(2LL*b - a)){
            long long R = (t + a)/2;
            vector<int> xs; xs.reserve((size_t)a);
            for(int i=1;i<=a;i++) xs.push_back(i);
            long long cur = 1LL*a*(a+1)/2;
            for(long long i=a-1;i>=0;i--){
                int maxv = b - (a-1 - i);
                long long add = min<long long>(maxv - xs[i], R - cur);
                if(add>0){ xs[i] += (int)add; cur += add; }
            }
            for(int i=n;i>b;i--) order.push_back(i);
            vector<char> chosen(b+1,false);
            for(int v: xs) chosen[v]=true;
            for(int v: xs) order.push_back(v);
            for(int i=b;i>=1;i--) if(!chosen[i]) order.push_back(i);

            for(size_t i=0;i<order.size();++i){
                if(i) cout << ' ';
                cout << (2LL*order[i]-1);
            }
            cout << '\n';
            return 0;
        }
    }

    cout << "impossible\n";
    return 0;
}