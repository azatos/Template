#include <bits/stdc++.h>
#define x first
#define y second
using namespace std;

typedef pair<double,double> P; // Point
struct Q{double x, y, r;}; // Cirlce

const double eps = 1e-9;

// find two point of regular triangle with given two point
vector<P> R_Tri(P a, P b){
}

// return distance between two points
double Dis(P a, P b){
	return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

// return intersections of two line
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
	
	cout << x1 << " " << y1 << " " << dx << " " << dy << " " << dr << endl;
	
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

int main(){
	int n;
	scanf("%d",&n);
	vector<P> v(n);
	for(auto &t:v) scanf("%lf %lf",&t.x,&t.y);
	
	auto p = CC({0,0,1},{2,0,1});
	cout << p.size() << endl;
	for(auto &t:p){
		cout << t.x << " " << t.y << endl;
	}
}

