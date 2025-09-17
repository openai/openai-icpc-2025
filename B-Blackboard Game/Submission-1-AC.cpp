#include <bits/stdc++.h>
using namespace std;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int t;
    if(!(cin>>t)) return 0;
    vector<int> ns(t);
    int maxn=0;
    for(int i=0;i<t;++i){ cin>>ns[i]; maxn=max(maxn, ns[i]); }

    const int SM=83;
    vector<int> small(SM+1);
    vector<int> arr = {0,0,0,2,0,0,0,0,2,0,0,0,2,2,2,0,4,4,0,0,2,0,0,0,4,4,4,4,0,0,2,2,2,2,2,0,4,4,4,4,4,4,4,4,4,0,0,0,2,2,2,0,2,2,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,0,4,0,0,2,0,0,0};
    for(int i=0;i<=SM && i<(int)arr.size();++i) small[i]=arr[i];

    int limit = max(0, maxn/3);
    vector<int> prevPrime(limit+1,0);
    if(maxn>=84){
        vector<char> is(limit+1,true);
        if(limit>=0) is[0]=false;
        if(limit>=1) is[1]=false;
        for(long long i=2;i*i<=limit;i++) if(is[i]) for(long long j=i*i;j<=limit;j+=i) is[j]=false;
        int last=0; 
        for(int i=2;i<=limit;i++){ if(is[i]) last=i; prevPrime[i]=last; }
    }

    for(int n: ns){
        if(n<=SM){
            int s = small[n];
            if(s==0) cout<<"second\n"; 
            else cout<<"first "<<s<<"\n";
        } else {
            int p = prevPrime[n/3];
            if(p==0) p=2;
            cout<<"first "<<2*p<<"\n";
        }
    }
    return 0;
}