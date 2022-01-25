#include <bits/stdc++.h>
#define x first
#define y second
using namespace std;

typedef pair<double,double> P; // Point
struct Q{double x, y, r;}; // Cirlce

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
			B = 2*(dx*(a.x-c.x) + dy*(a.y-c.y)),
			C = (a.x-c.x)*(a.x-c.x) + (a.y-c.y)*(a.y-c.y) - c.r*c.r;
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
vector<P> CC(Q a, Q b){
	double dx = 2*(a.y-b.y), dy = -2*(a.x-b.x), dr = (a.x*a.x-b.x*b.x + a.y*a.y-b.y*b.y - a.r*a.r+b.r*b.r);
	double x1 = (dx==0) ? -dr/dy : 0, y1 = (dx==0) ? 0 : (dr+dy*x1)/dx;
	double	A = dx*dx + dy*dy,
			B = 2*(dx*(x1-a.x) + dy*(y1-a.y)),
			C = (x1-a.x)*(x1-a.x) + (y1-a.y)*(y1-a.y) - a.r*a.r;
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

// return convex hull with graham scan
vector<P> ConvexHull(vector<P> v){
	if(v.size()<3) return {};
	sort(v.begin(), v.end());
	sort(v.begin()+1, v.end(), [&](P x,P y){ return ccw(v[0],x,y) > 0 || (ccw(v[0],x,y) == 0 && Dis(v[0],x) < Dis(v[0],y)); });
	vector<P> r{v[0],v[1]};
	for(int i=2; i<v.size(); r.push_back(v[i]), ++i) while(r.size()>1 && ccw(r[r.size()-2],r[r.size()-1],v[i])<=0) r.pop_back();
	return r;
}

// return MST length of ConvexHull
double MST_D(vector<P> v){
	double D=Dis(v[0],v[v.size()-1]), M = D;
	for(int i=1; i<v.size(); ++i) D+=Dis(v[i-1],v[i]), M=max(M,Dis(v[i-1],v[i]));
	return D-M;
}

// Triangle Steiner Points
vector<P> TSP(vector<P> v){
	vector<P> t;
	t = CC({v[0].x,v[0].y,Dis(v[0],v[1])},{v[1].x,v[1].y,Dis(v[0],v[1])});
	P X = ccw(v[0],v[1],v[2]) == ccw(v[0],v[1],t[0]) ? t[1] : t[0];
	t = CC({v[1].x,v[1].y,Dis(v[1],v[2])},{v[2].x,v[2].y,Dis(v[1],v[2])});
	P Y = ccw(v[1],v[2],v[0]) == ccw(v[1],v[2],t[0]) ? t[1] : t[0];
	return LL(X,v[2],Y,v[0]);
}
// Length of Triangle Steiner Tree
double TSP_D(vector<P> v){
	P Z = TSP(v)[0];
	return Dis(Z,v[0]) + Dis(Z,v[1]) + Dis(Z,v[2]);
}

vector<P> QSP(vector<P> v){
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
		puts("TSP");
		for(auto &t:u) printf("%lf %lf\n",t.x,t.y);
		printf("TSP : %lf\t MST : %lf\n",TSP_D(w),MST_D(w));
	}
}

