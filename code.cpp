#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

using namespace std;

class term {
public:
    long long a;
    int b, c, d;

    term() : a(0), b(0), c(0), d(0) {}
    term(long long _a, int _b, int _c, int _d) : a(_a), b(_b), c(_c), d(_d) {}

    bool operator==(const term& obj) const {
        return b == obj.b && c == obj.c && d == obj.d;
    }
    bool operator!=(const term& obj) const {
        return !(*this == obj);
    }
    bool operator<(const term& other) const {
        if (b != other.b) return b > other.b;
        if (c != other.c) return c > other.c;
        return d > other.d;
    }
};

class poly {
public:
    int n;
    term *t;

    poly() : n(0), t(nullptr) {}
    poly(int _n) : n(_n) {
        t = n > 0 ? new term[n] : nullptr;
    }
    poly(const poly &p) : n(p.n) {
        t = n > 0 ? new term[n] : nullptr;
        for (int i = 0; i < n; ++i) t[i] = p.t[i];
    }
    poly(long long a) {
        if (a != 0) {
            n = 1;
            t = new term[1];
            t[0] = term(a, 0, 0, 0);
        } else {
            n = 0;
            t = nullptr;
        }
    }
    poly(const term& _t) {
        if (_t.a != 0) {
            n = 1;
            t = new term[1];
            t[0] = _t;
        } else {
            n = 0;
            t = nullptr;
        }
    }

    void simplify() {
        if (n == 0) return;
        sort(t, t + n);
        int new_n = 0;
        for (int i = 0; i < n; ++i) {
            if (t[i].a == 0) continue;
            if (new_n > 0 && t[i] == t[new_n - 1]) {
                t[new_n - 1].a += t[i].a;
                if (t[new_n - 1].a == 0) new_n--;
            } else {
                t[new_n++] = t[i];
            }
        }
        n = new_n;
    }

    poly operator+(const poly &obj) const {
        poly ans(n + obj.n);
        for (int i = 0; i < n; ++i) ans.t[i] = t[i];
        for (int i = 0; i < obj.n; ++i) ans.t[i + n] = obj.t[i];
        ans.simplify();
        return ans;
    }

    poly operator-(const poly &obj) const {
        poly ans(n + obj.n);
        for (int i = 0; i < n; ++i) ans.t[i] = t[i];
        for (int i = 0; i < obj.n; ++i) {
            ans.t[i + n] = obj.t[i];
            ans.t[i + n].a *= -1;
        }
        ans.simplify();
        return ans;
    }

    poly operator*(const poly &obj) const {
        if (n == 0 || obj.n == 0) return poly(0LL);
        poly ans(n * obj.n);
        int k = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < obj.n; ++j) {
                ans.t[k++] = term(t[i].a * obj.t[j].a, t[i].b + obj.t[j].b, t[i].c + obj.t[j].c, t[i].d + obj.t[j].d);
            }
        }
        ans.n = k;
        ans.simplify();
        return ans;
    }

    poly& operator=(const poly &obj) {
        if (&obj == this) return *this;
        delete[] t;
        n = obj.n;
        t = n > 0 ? new term[n] : nullptr;
        for (int i = 0; i < n; ++i) t[i] = obj.t[i];
        return *this;
    }

    poly derivate() const {
        poly ans(3 * n);
        int k = 0;
        for (int i = 0; i < n; ++i) {
            if (t[i].b > 0) ans.t[k++] = term(t[i].a * t[i].b, t[i].b - 1, t[i].c, t[i].d);
            if (t[i].c > 0) ans.t[k++] = term(t[i].a * t[i].c, t[i].b, t[i].c - 1, t[i].d + 1);
            if (t[i].d > 0) ans.t[k++] = term(-t[i].a * t[i].d, t[i].b, t[i].c + 1, t[i].d - 1);
        }
        ans.n = k;
        ans.simplify();
        return ans;
    }

    string to_string_poly() const {
        if (n == 0) return "0";
        string res = "";
        for (int i = 0; i < n; ++i) {
            long long a = t[i].a;
            if (a > 0 && i > 0) res += "+";
            if (a < 0) {
                res += "-";
                a = -a;
            }
            bool is_const = (t[i].b == 0 && t[i].c == 0 && t[i].d == 0);
            if (a != 1 || is_const) {
                res += to_string(a);
            }
            if (t[i].b > 0) {
                res += "x";
                if (t[i].b > 1) res += "^" + to_string(t[i].b);
            }
            if (t[i].c > 0) {
                res += "sin";
                if (t[i].c > 1) res += "^" + to_string(t[i].c);
                res += "x";
            }
            if (t[i].d > 0) {
                res += "cos";
                if (t[i].d > 1) res += "^" + to_string(t[i].d);
                res += "x";
            }
        }
        return res;
    }

    ~poly() {
        delete[] t;
    }
};

class frac {
public:
    poly p, q;

    frac() {}
    frac(long long val) : p(val), q(1LL) {}
    frac(const term& _t) : p(_t), q(1LL) {}
    frac(const poly& _p, const poly& _q) : p(_p), q(_q) {}

    frac operator+(const frac &obj) const {
        return frac(p * obj.q + q * obj.p, q * obj.q);
    }
    frac operator-(const frac &obj) const {
        return frac(p * obj.q - q * obj.p, q * obj.q);
    }
    frac operator*(const frac &obj) const {
        return frac(p * obj.p, q * obj.q);
    }
    frac operator/(const frac &obj) const {
        return frac(p * obj.q, q * obj.p);
    }
    frac derivate() const {
        return frac(p.derivate() * q - q.derivate() * p, q * q);
    }
    void output() const {
        if (p.n == 0) {
            cout << "0" << endl;
            return;
        }
        bool q_is_one = (q.n == 1 && q.t[0].a == 1 && q.t[0].b == 0 && q.t[0].c == 0 && q.t[0].d == 0);
        if (q_is_one) {
            cout << p.to_string_poly() << endl;
            return;
        }
        string sp = p.to_string_poly();
        string sq = q.to_string_poly();
        if (p.n > 1) cout << "(" << sp << ")";
        else cout << sp;
        cout << "/";
        if (q.n > 1) cout << "(" << sq << ")";
        else cout << sq;
        cout << endl;
    }
};

class Parser {
    string s;
    int pos;

public:
    Parser(string _s) : s(_s), pos(0) {}

    frac parseExpression() {
        frac res = parseTerm();
        while (pos < (int)s.length() && (s[pos] == '+' || s[pos] == '-')) {
            char op = s[pos++];
            frac next = parseTerm();
            if (op == '+') res = res + next;
            else res = res - next;
        }
        return res;
    }

    frac parseTerm() {
        frac res = parseFactor();
        while (pos < (int)s.length() && (s[pos] == '*' || s[pos] == '/')) {
            char op = s[pos++];
            frac next = parseFactor();
            if (op == '*') res = res * next;
            else res = res / next;
        }
        return res;
    }

    frac parseFactor() {
        if (pos < (int)s.length() && s[pos] == '(') {
            pos++;
            frac res = parseExpression();
            if (pos < (int)s.length() && s[pos] == ')') pos++;
            return res;
        } else if (pos < (int)s.length() && s[pos] == '-') {
            pos++;
            return frac(-1LL) * parseFactor();
        } else {
            return parseBaseTerm();
        }
    }

    frac parseBaseTerm() {
        long long a = 1;
        bool has_a = false;
        if (pos < (int)s.length() && isdigit(s[pos])) {
            a = 0;
            while (pos < (int)s.length() && isdigit(s[pos])) {
                a = a * 10 + (s[pos++] - '0');
            }
            has_a = true;
        }

        int b = 0, c = 0, d = 0;
        while (pos < (int)s.length() && (s[pos] == 'x' || s[pos] == 's' || s[pos] == 'c')) {
            if (s[pos] == 'x') {
                pos++;
                int exp = 1;
                if (pos < (int)s.length() && s[pos] == '^') {
                    pos++;
                    exp = 0;
                    while (pos < (int)s.length() && isdigit(s[pos])) {
                        exp = exp * 10 + (s[pos++] - '0');
                    }
                }
                b += exp;
            } else if (pos + 3 <= (int)s.length() && s.substr(pos, 3) == "sin") {
                pos += 3;
                int exp = 1;
                if (pos < (int)s.length() && s[pos] == '^') {
                    pos++;
                    exp = 0;
                    while (pos < (int)s.length() && isdigit(s[pos])) {
                        exp = exp * 10 + (s[pos++] - '0');
                    }
                }
                if (pos < (int)s.length() && s[pos] == 'x') pos++;
                c += exp;
            } else if (pos + 3 <= (int)s.length() && s.substr(pos, 3) == "cos") {
                pos += 3;
                int exp = 1;
                if (pos < (int)s.length() && s[pos] == '^') {
                    pos++;
                    exp = 0;
                    while (pos < (int)s.length() && isdigit(s[pos])) {
                        exp = exp * 10 + (s[pos++] - '0');
                    }
                }
                if (pos < (int)s.length() && s[pos] == 'x') pos++;
                d += exp;
            } else {
                break;
            }
        }
        return frac(term(a, b, c, d));
    }
};

void solve(char *s, int n) {
    Parser p(s);
    frac f = p.parseExpression();
    f.output();
    frac df = f.derivate();
    df.output();
}

int main() {
    string str;
    if (!(cin >> str)) return 0;
    int n = str.length();
    char *s = new char[n + 2]{0};
    for (int i = 0; i < n; ++i) s[i] = str[i];
    solve(s, n);
    delete []s;
    return 0;
}
