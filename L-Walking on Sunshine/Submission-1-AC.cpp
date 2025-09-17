#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int n;
    long long xc, yc, xa, ya;
    if(!(cin >> n >> xc >> yc >> xa >> ya)) return 0;

    if (ya >= yc) {
        cout.setf(std::ios::fixed);
        cout << setprecision(10) << 0.0 << "\n";
        // still need to read remaining rectangles to avoid premature EOF?
        for(int i=0;i<n;i++){
            long long x1,y1,x2,y2; cin>>x1>>y1>>x2>>y2;
        }
        return 0;
    }
    
    long long ylow = ya, yhigh = yc;
    vector<pair<long long,long long>> segs;
    segs.reserve(n);
    for (int i = 0; i < n; ++i) {
        long long x1,y1,x2,y2;
        cin >> x1 >> y1 >> x2 >> y2;
        long long l = max(ylow, y1);
        long long r = min(yhigh, y2);
        if (l < r) segs.push_back({l, r});
    }

    sort(segs.begin(), segs.end());
    long long unionLen = 0;
    bool has = false;
    long long curL = 0, curR = -1;
    for (auto &p : segs) {
        long long l = p.first, r = p.second;
        if (!has) { curL = l; curR = r; has = true; }
        else if (l <= curR) curR = max(curR, r);
        else { unionLen += curR - curL; curL = l; curR = r; }
    }
    if (has) unionLen += curR - curL;

    long double res = (long double)(yhigh - ylow) - (long double)unionLen;
    if (res < 0) res = 0;
    cout.setf(std::ios::fixed);
    cout << setprecision(10) << (double)res << "\n";
    return 0;
}