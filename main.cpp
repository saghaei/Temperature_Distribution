/* Temperature distribution in the rectangular plate */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

void print(double**, const size_t&, const size_t&);
bool test(double**, double**, const size_t&, const size_t&);
void pointGaussSeidel(double**, double**, const size_t&, const size_t&, const double&);
void lineGaussSeidel(double**, double**, const size_t&, const size_t&, const double&);
void alternatingDirectImplicit(double**, double**, const size_t&, const size_t&, const double&, double&);
void setWeight(double&);
size_t counter=0;

int main() {
double u1=100, u2=200, u3=300, u4=400;
size_t r_max=5, c_max=5;
double b=1,w=1;
unsigned short int mtd;

cout << "Steady-State temperature distribution on a two-dimensional\n";
cout << "rectangular plate with dirichlet boundary condition.\n" << endl;
cout << "Enter Up/Left/Down/Right temperature consequently: ";
cin >> u1 >> u2 >> u3 >> u4;
cout << "Enter Row/Column consequently: ";
cin >> r_max >> c_max;

double **u = new double*[r_max];
for(size_t c=0 ; c<r_max ; c++)
u[c] = new double[c_max];
double **p = new double*[r_max];
for(size_t c=0 ; c<r_max ; c++)
p[c] = new double[c_max];

for(size_t c=1 ; c<c_max-1 ; c++)
u[0][c] = u1;
for(size_t r=1 ; r<r_max-1 ; r++)
u[r][0] = u2;
for(size_t c=1 ; c<c_max-1 ; c++)
u[r_max-1][c] = u3;
for(size_t r=1 ; r<r_max-1 ; r++)
u[r][c_max-1] = u4;

for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
u[r][c] = 0;
for(size_t r=0 ; r<r_max ; r++)
for(size_t c=0 ; c<c_max ; c++)
p[r][c] = u[r][c];

cout << "1) Point Gauss-Seidel\n2) Line Gauss-Seidel\n3) Alternative Direct Implicit\n" << endl;

do {
cout << "Enter solving method: ";
cin >> mtd;
}while(mtd<1 || mtd>5);

cout << endl;

switch (mtd) {
case 1:
pointGaussSeidel(u,p,r_max,c_max,b);
break;
case 2:
lineGaussSeidel(u,p,r_max,c_max,b);
break;
case 3:
alternatingDirectImplicit(u,p,r_max,c_max,b,w);
break;
default:
cout << "Program execution ended." << endl;
exit(0);
}

print(u, r_max, c_max);
fstream output("data", ios::out);
if(!output) {
cerr << "File could not be created." << endl;
exit(1);
}
for(size_t r=0 ; r<r_max ; r++) {
for(size_t c=0 ; c<c_max ; c++)
if((r==0 && c==0) || (r==0 && c==c_max-1) || (r==r_max-1 && c==0) || (r==r_max-1 && c==c_max-1))
output << setw(7) << left << '*' << '\t';
else
output << setw(7) << setprecision(2) << fixed << left << u[r][c] << '\t';
output << endl;
}

for(size_t r=0 ; r<r_max ; r++)
delete[] p[r];
delete[] p;
for(size_t r=0 ; r<r_max ; r++)
delete[] u[r];
delete[] u;
cout << endl;
return 0;
}

void print(double** u, const size_t& r_max, const size_t& c_max) {
for(size_t r=0 ; r<r_max ; r++) {
for(size_t c=0 ; c<c_max ; c++)
if((r==0 && c==0) || (r==0 && c==c_max-1) || (r==r_max-1 && c==0) || (r==r_max-1 && c==c_max-1))
cout << setw(6) << left << '*' << '\t';
else
cout << setprecision(1) << fixed << setw(6) << u[r][c] << '\t';
cout << endl;
}
}

bool test(double** u, double** p, const size_t& r_max, const size_t& c_max) {
for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
if(fabs(u[r][c] - p[r][c]) > 0.01)
return true;
return false;
}

void pointGaussSeidel(double** u, double** p, const size_t& r_max, const size_t& c_max, const double& b) {
do {
for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
p[r][c] = u[r][c];
for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
u[r][c] = (0.5/(1+b*b)) * (u[r+1][c] + u[r-1][c] + b*b*(u[r][c+1] + u[r][c-1]));
counter++;
} while (test(u,p,r_max,c_max));
cout << "Iterations: " << counter << endl << endl;
counter=0;
}

void lineGaussSeidel(double** u, double** p, const size_t& r_max, const size_t& c_max, const double& b) {
do {
for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
p[r][c] = u[r][c];
for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
u[r][c] = 0.5*(u[r-1][c] + u[r+1][c] + b*b*p[r][c+1] + b*b*u[r][c-1])/(1+b*b);
counter++;
} while (test(u,p,r_max,c_max));
cout << "Iterations: " << counter << endl << endl;
counter=0;
}

void alternatingDirectImplicit(double** u, double** p, const size_t& r_max, const size_t& c_max, const double& b, double& w) {
setWeight(w);
double **m = new double*[r_max];
for(size_t c=0 ; c<r_max ; c++)
m[c] = new double[c_max];
for(size_t r=0 ; r<r_max ; r++)
for(size_t c=0 ; c<c_max ; c++)
m[r][c] = u[r][c];
do {
for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
p[r][c] = u[r][c];
for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
m[r][c] = (0.5/(1+b*b))*(w*(m[r-1][c]+m[r+1][c])+2*(1-w)*(1+b*b)*p[r][c]+w*b*b*(p[r][c+1]+m[r][c-1]));
for(size_t r=1 ; r<r_max-1 ; r++)
for(size_t c=1 ; c<c_max-1 ; c++)
u[r][c] = (0.5/(1+b*b))*(w*b*b*(u[r][c-1]+u[r][c+1])+2*(1-w)*(1+b*b)*m[r][c]+w*(m[r+1][c]+u[r-1][c]));
counter++;
} while(test(u,p,r_max,c_max));
for(size_t r=0 ; r<r_max ; r++)
delete[] m[r];
delete[] m;
cout << "Iterations: " << counter << endl << endl;
counter=0;
}

void setWeight(double &w) {
do {
cout << "Enter (0<w<2): ";
cin >> w;
} while(w<=0 || w>=2);
}
