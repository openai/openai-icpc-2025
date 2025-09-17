#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int r, c;
    if (!(cin >> r >> c)) return 0;
    vector<string> grid(r);
    for (int i = 0; i < r; ++i) cin >> grid[i];
    string s;
    cin >> s;
    int n = s.size();

    pair<int,int> start = {-1, -1};
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            if (grid[i][j] == 'S') start = {i, j};

    vector<int> diridx(256, -1);
    diridx['N'] = 0; diridx['E'] = 1; diridx['S'] = 2; diridx['W'] = 3;
    int dr[4] = {-1, 0, 1, 0};
    int dc[4] = {0, 1, 0, -1};

    vector<array<int,4>> perms;
    array<int,4> base = {0, 1, 2, 3};
    sort(base.begin(), base.end());
    do { perms.push_back(base); } while (next_permutation(base.begin(), base.end()));
    int P = perms.size();
    vector<array<int,4>> pos(P);
    for (int i = 0; i < P; ++i)
        for (int k = 0; k < 4; ++k)
            pos[i][perms[i][k]] = k;

    vector<vector<char>> allowed(n, vector<char>(P, false));

    int r0 = start.first, c0 = start.second;
    for (int i = 0; i < n; ++i) {
        bool valid[4];
        for (int d = 0; d < 4; ++d) {
            int nr = r0 + dr[d], nc = c0 + dc[d];
            if (nr < 0 || nr >= r || nc < 0 || nc >= c) valid[d] = false;
            else valid[d] = (grid[nr][nc] != '#');
        }
        int sd = diridx[(int)s[i]];
        for (int p = 0; p < P; ++p) {
            bool ok = true;
            for (int d = 0; d < 4; ++d) if (d != sd && valid[d]) {
                if (pos[p][sd] >= pos[p][d]) { ok = false; break; }
            }
            allowed[i][p] = ok;
        }
        r0 += dr[sd]; c0 += dc[sd];
    }

    const int INF = 1e9;
    vector<int> dpPrev(P, INF), dpCur(P, INF);
    for (int p = 0; p < P; ++p) dpPrev[p] = allowed[0][p] ? 0 : INF;

    for (int i = 1; i < n; ++i) {
        fill(dpCur.begin(), dpCur.end(), INF);
        for (int q = 0; q < P; ++q) if (dpPrev[q] < INF) {
            for (int p = 0; p < P; ++p) if (allowed[i][p]) {
                int cost = dpPrev[q] + (p != q);
                if (cost < dpCur[p]) dpCur[p] = cost;
            }
        }
        dpPrev.swap(dpCur);
    }

    int ans = *min_element(dpPrev.begin(), dpPrev.end());
    cout << ans << "\n";
    return 0;
}