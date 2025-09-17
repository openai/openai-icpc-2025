#include <bits/stdc++.h>
using namespace std;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; if(!(cin>>n)) return 0;
    int k; cin>>k;

    vector<int> pos(n+1,0);

    auto rotate = [&](int i, int delta){
        cout<<i<<" "<<delta<<"\n"<<flush;
        pos[i] = ((pos[i] + delta) % n + n) % n;
        if(!(cin>>k)) exit(0);
        if(k==1) exit(0);
        return k;
    };

    auto rotateTo = [&](int i, int target){
        int cur = pos[i];
        int diff = target - cur;
        diff %= n; if(diff < 0) diff += n;
        if(diff > n - diff) diff -= n;
        if(diff != 0) rotate(i, diff);
        return k;
    };

    auto scan = [&](int i){
        vector<int> arr(n);
        arr[0]=k;
        for(int t=1;t<=n-1;++t){
            rotate(i,1);
            arr[t]=k;
        }
        rotate(i,1);
        int c=*min_element(arr.begin(), arr.end());
        vector<int> mem(n);
        for(int t=0;t<n;++t) mem[t] = (arr[t]==c) ? 1 : 0;
        return mem;
    };

    vector<int> mem1 = scan(1);
    int s_index=-1, r_index=-1;
    for(int t=0;t<n;++t){
        if(mem1[t]==0 && s_index==-1) s_index=t;
        if(mem1[t]==1 && r_index==-1) r_index=t;
    }
    if(s_index==-1) s_index=0;
    if(r_index==-1) r_index=(s_index+1)%n;

    rotateTo(1, s_index);
    vector<int> spos(n+1,-1);
    spos[1]=s_index;

    for(int j=2;j<=n;++j){
        vector<int> Awith = scan(j);
        rotateTo(1, r_index);
        vector<int> Awithout = scan(j);
        int found=-1;
        for(int t=0;t<n;++t){
            if(Awith[t]==1 && Awithout[t]==0){ found=t; break; }
        }
        if(found==-1){
            vector<int> diffIdx;
            for(int t=0;t<n;++t) if(Awith[t]!=Awithout[t]) diffIdx.push_back(t);
            if(!diffIdx.empty()) found=diffIdx[0]; else found=0;
        }
        spos[j]=found;
        rotateTo(1, s_index);
    }

    for(int i=1;i<=n;++i){
        rotateTo(i, spos[i]);
        if(k==1) return 0;
    }
    return 0;
}