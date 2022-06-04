
#include "iostream"
#include "cmath"
#include "iomanip"
#include "fstream"

using namespace std;

int main(){
    // deklarasi variabel
    float pi = 3.14;
    
    // deklarasi molekul air
    float Mgx, Mgy, Mgz;
    float Cl1x, Cl1y, Cl1z;
    float Cl2x, Cl2y, Cl2z;

    // masukkan nilai untuk setiap variabel MgCl2
    Mgx = Mgy = Mgz = 0.0;
    Cl1x = sin(88.75*pi/180) * 177.5;
    Cl1y = 0.0;
    Cl1z = cos(88.75*pi/180) * 177.5;

    Cl2x = - Cl1x;
    Cl2y = Cl1y;
    Cl2z = Cl1z;

    // memasukkan nilai tersebut ke dalam file
    ofstream file;
    file.open("MgCl2.xyz")
    file << "3\n" << endl;
    file << setw(3) << "Mg" << setw(3) << " " 
    << fixed << setprecision(3) << Mgx << setw(3) << " " \
    << fixed << setprecision(3) << Mgy << setw(3) << " " \
    << fixed << setprecision(3) << Mgz << "\n";
