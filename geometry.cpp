#include <bits/stdc++.h>
#define x first
#define y second
using namespace std;

typedef pair<double,double> P; // Point
struct Q{P p; double r;}; // Cirlce

vector<P> mergeP(vector<P> a, vector<P> b){
	a.insert(a.end(),b.begin(),b.end());
	return a;
}

const double eps = 1e-9;
int ccw(P a, P b, P c){
	double X = a.x*b.y + b.x*c.y + c.x*a.y;
	double Y = a.y*b.x + b.y*c.x + c.y*a.x;
	return X>Y?1:X<Y?-1:0;
}

// return distance between two points
double Dis(P a, P b){
	return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

// return intersections of two line _ ab | cd
vector<P> LL(P a, P b, P c, P d){
	if((a.y - b.y)*(c.x - d.x) == (c.y - d.y)*(a.x - b.x)) return {};
	return {{((a.x*b.y - a.y*b.x)*(c.x-d.x) - (a.x-b.x)*(c.x*d.y - c.y*d.x)) / ((a.x - b.x)*(c.y - d.y) - (a.y - b.y)*(c.x - d.x)),
			((a.x*b.y - a.y*b.x)*(c.y-d.y) - (a.y-b.y)*(c.x*d.y - c.y*d.x)) / ((a.x - b.x)*(c.y - d.y) - (a.y - b.y)*(c.x - d.x))}};
}

// return intersections of circle and line
vector<P> LC(P a, P b, Q c){
	double dx = (a.x-b.x), dy = (a.y-b.y);
	double	A = dx*dx + dy*dy,
			B = 2*(dx*(a.x-c.p.x) + dy*(a.y-c.p.y)),
			C = (a.x-c.p.x)*(a.x-c.p.x) + (a.y-c.p.y)*(a.y-c.p.y) - c.r*c.r;
	double	D = B*B - 4*A*C;
	if(D<0 || A<=eps) return {};
	if(D==0){
		double T = -B/(2*A);
		return {{a.x+T*dx, a.y+T*dy}};
	}
	double T1 = (-B - sqrt(D))/(2*A), T2 = (-B + sqrt(D))/(2*A);
	return {{a.x+T1*dx, a.y+T1*dy},
			{a.x+T2*dx, a.y+T2*dy}};
}

// return intersections of two circle
vector<P> CC(Q X, Q Y){
	P a = X.p, b = Y.p;
	double dx = 2*(a.y-b.y), dy = -2*(a.x-b.x), dr = (a.x*a.x-b.x*b.x + a.y*a.y-b.y*b.y - X.r*X.r+Y.r*Y.r);
	double x1 = (dx==0) ? -dr/dy : 0, y1 = (dx==0) ? 0 : (dr+dy*x1)/dx;
	double	A = dx*dx + dy*dy,
			B = 2*(dx*(x1-a.x) + dy*(y1-a.y)),
			C = (x1-a.x)*(x1-a.x) + (y1-a.y)*(y1-a.y) - X.r*X.r;
	double	D = B*B - 4*A*C;
	if(D<0 || A<=eps) return {};
	if(D==0){
		double T = -B/(2*A);
		return {{x1+T*dx, y1+T*dy}};
	}
	double T1 = (-B - sqrt(D))/(2*A), T2 = (-B + sqrt(D))/(2*A);
	return {{x1+T1*dx, y1+T1*dy},
			{x1+T2*dx, y1+T2*dy}};
}

// return circumscribed circle of three point
Q OC(P a,P b,P c){
	vector<P> t1 = CC({a,Dis(a,b)},{b,Dis(a,b)}), t2 = CC({a,Dis(a,c)},{c,Dis(a,c)}), Z = LL(t1[0],t1[1],t2[0],t2[1]);
	return {Z[0],Dis(Z[0],a)};
}

// return convex hull with graham scan
vector<P> ConvexHull(vector<P> v){
	if(v.size()<3) return {};
	sort(v.begin(), v.end());
	sort(v.begin()+1, v.end(), [&](P x,P y){ return ccw(v[0],x,y) > 0 || (ccw(v[0],x,y) == 0 && Dis(v[0],x) < Dis(v[0],y)); });
	vector<P> r{v[0],v[1]};
	for(int i=2; i<v.size(); r.push_back(v[i]), ++i) while(r.size()>1 && ccw(r[r.size()-2],r[r.size()-1],v[i])<=0) r.pop_back();
	return r;
}

double MST(vector<P> v){
	vector<bool> vis(v.size(),0);
	vector<pair<double,pair<int,int>>> e;
	for(int i=0; i<v.size(); ++i) for(int j=i+1; j<v.size(); ++j) e.push_back({Dis(v[i],v[j]),{i,j}});
	sort(e.begin(),e.end());
	double ret = 0;
	for(auto [a,b]:e){
		if(vis[b.x] && vis[b.y]) continue;
		ret += a;
		vis[b.x] = vis[b.y] = 1;
	}
	return ret;
}

// return MST length of ConvexHull
double MST_D(vector<P> v){
	double D=Dis(v[0],v[v.size()-1]), M = D;
	for(int i=1; i<v.size(); ++i) D+=Dis(v[i-1],v[i]), M=max(M,Dis(v[i-1],v[i]));
	return D-M;
}

// Triangular Steiner Points
vector<vector<P>> TSP(vector<P> v){
	vector<P> t;
	t = CC({v[0],Dis(v[0],v[1])},{v[1],Dis(v[0],v[1])});
	P X = ccw(v[0],v[1],v[2]) == ccw(v[0],v[1],t[0]) ? t[1] : t[0];
	t = CC({v[1],Dis(v[1],v[2])},{v[2],Dis(v[1],v[2])});
	P Y = ccw(v[1],v[2],v[0]) == ccw(v[1],v[2],t[0]) ? t[1] : t[0];
	return {LL(X,v[2],Y,v[0])};
}

// Quadrilateral Steiner Points
vector<vector<P>> QSP(vector<P> v){
	vector<vector<P>> r;
	vector<tuple<int,int,int,int>> T = {{0,1,2,3},{1,2,3,0}};
	for(auto &[a,b,c,d]:T){
		vector<P> t;
		P X, Y, A, B;
		t = CC({v[a],Dis(v[a],v[b])},{v[b],Dis(v[a],v[b])});
		X = ccw(v[a],v[b],v[c]) == ccw(v[a],v[b],t[0]) ? t[1] : t[0];
		t = CC({v[c],Dis(v[c],v[d])},{v[d],Dis(v[c],v[d])});
		Y = ccw(v[c],v[d],v[a]) == ccw(v[c],v[d],t[0]) ? t[1] : t[0];
				
		t = LC(X,Y,OC(v[a],v[b],X));
		A = Dis(X,t[0]) < Dis(X,t[1]) ? t[1] : t[0];
		t = LC(X,Y,OC(v[c],v[d],Y));
		B = Dis(Y,t[0]) < Dis(Y,t[1]) ? t[1] : t[0];
		r.push_back({A,B});
	}
	return r;
}

// Pentagonal Steiner Points
vector<vector<P>> PSP(vector<P> v){
}

int main(){
	int n;
	scanf("%d",&n);
	vector<P> v(n);
	for(auto &t:v) scanf("%lf %lf",&t.x,&t.y);
	auto w = ConvexHull(v);
	if(w.size()!=n) { puts("Error! Concave!"); return 0; }
	
	if(n==3){
		auto u = TSP(w);
		printf("MST : %lf\n",MST(w));
		for(auto&tt:u) {
			for(auto &t:tt) printf("%lf %lf\n",t.x,t.y);
			printf("TSP : %lf\n",MST(mergeP(w,tt)));
		}
	}
	else if(n==4){
		auto u = QSP(w);
		printf("MST : %lf\n",MST(w));
		for(auto&tt:u) {
			for(auto &t:tt) printf("%lf %lf\n",t.x,t.y);
			printf("QSP : %lf\n",MST(mergeP(w,tt)));
		}
	}
	else if(n==5){
	}
}

