gws rahma 
#sertakan "iostream
#sertakan "cmath"
#include "iomanip"
#include "fstream"

using namespace std;

float densitas;         /*densitas ini memiliki satuan g cm^-3
                          dan ini harus diinput dari user */
// converter dari nilai cm^-3 ke A^3
float cm3_to_A3 = 1.0E-24;
float numb_mol;
float Na = 6.022E23;            // bilangan avogadro
float Mr_metana = 16.0;         // massa relatif dari metana
int numb_lat;                   // panjang sel simulasi
int N;                          // jumlah molekul
float volum;                    // volume dari molekul metana
float lx, ly, lz;               // panjang sel simulasi terhitung 
float lat;                      // jarak antar molekul metana

// variabel perhitungan energi potensial
float sigma = 2.725;            // satuan (A)
float epsilon = 3.725;          // satuan (kj/mol)
float rij, rcut, rcut2;
float Ep, Ep_LJ, A12, B6;
float rij6, rij12;
float dx, dy, dz;

int main(){
    
    rcut = 3 * sigma;
    rcut2 = pow(rcut,2);
    
    cout << "program menghitung energi potensial dari metana" << endl;
    cout << "Masukkan nilai densitas dari metana (eg. 1.0): ";
    cin >> densitas;

    numb_mol = densitas * (Na / Mr_metana) * cm3_to_A3;
    cout << "Masukkan panjang sel simulasi (eg. 5): ";
    cin >> numb_lat;
    N = (pow(numb_lat,3)) * 3;
    volum = (float)N / numb_mol;

    // panjang sel simulasi secara perhitungan
    lx = pow(volum,(1.0/3.0));
    ly = lx;
    lz = lx;
    if(lx < (2*rcut)){
        cout << "sorry mazeh, sel simulasi kamu besar";
        exit(0);
    } else{
        cout << "panjang sek simulasi baru: " << lx << endl;
    }
    lat = lx / (float)numb_lat;

    // deklarasi variabel
    float pi = 3.14;
    
    // deklarasi molekul metana
    float Cx[125], Cy[125], Cz[125];
    float H1x[125], H1y[125], H1z[125];
    float H2x[125], H2y[125], H2z[125];
    float H3x[125], H3y[125], H3z[125];
    float H4x[125], H4y[125], H4z[125];

    /*    H
          |
          C
        / | \
      H   H   H
     molekul air itu bentuk rigid
   */

   float rB_Cx = 0.0;
   float rB_Cy = 0.0;
   float rB_Cz = 0.0;

   float rB_H1x = sin(109.5*pi/180.0) * 1.090;
   float rB_H1y = 0.0;
   float rB_H1z = cos(109.5*pi/180.0) * 1.090;

   float rB_H2x = - rB_H1x;
   float rB_H2y = rB_H1y;
   float rB_H2z = rB_H1z;

   float rB_H3x = 0.0;
   float rB_H3y = cos(54.75*pi/180) * 1.090;
   float rB_H3z = - sin(54.75*pi/180) * 1.090;

   float rB_H4x = rB_H3x;
   float rB_H4y = - rB_H3y;
   float rB_H4z = rB_H3z;

    int counter = 0;
        for(int i = 0; i < numb_lat; i++){ //x
           for(int j = 0; j < numb_lat; j++){ //y
               for(int k = 0; k < numb_lat; k++){ //z 
                
                Cx[counter] = rB_Cx + (i*lat);
                Cy[counter] = rB_Cy + (j*lat);
                Cz[counter] = rB_Cz + (k*lat);

                H1x[counter] = Cx[counter] + rB_H1x;
                H1y[counter] = Cy[counter] + rB_H1y;
                H1z[counter] = Cz[counter] + rB_H1z;

                H2x[counter] = Cx[counter] + rB_H2x;
                H2y[counter] = Cy[counter] + rB_H2y;
                H2z[counter] = Cz[counter] + rB_H2z;

                H3x[counter] = Cx[counter] + rB_H3x;
                H3y[counter] = Cy[counter] + rB_H3y;
                H3z[counter] = Cz[counter] + rB_H3z;

                H4x[counter] = Cx[counter] + rB_H4x;
                H4y[counter] = Cy[counter] + rB_H4y;
                H4z[counter] = Cz[counter] + rB_H4z;

                counter += 1;
           }
        }
    }
    // main hitung energi potensial
    A12 = 4.0 * epsilon * pow(sigma,12);
    B6 = 4.0 * epsilon * pow(sigma, 6);

    Ep = 0.0;
    for(int a = 0; a < (N-1); a++){
        for(int b = (0+1); b < N; b++){
            
            dx = Cx[a] - Cx[b];
            dy = Cy[a] - Cy[b];
            dz = Cz[a] - Cz[b];

            dx = dx - round(dx/lx) * lx;
            dy = dy - round(dy/ly) * ly;
            dz = dz - round(dz/lz) * lz;

            rij = pow(dx,2) + pow(dy,2) + pow(dz,2);
            if(rij < rcut2){
                rij6 = pow(rij,3);
                rij12 = pow(rij6,2);
                Ep = (A12/rij12) - (B6/rij6);
                Ep_LJ = (A12/rij12) - (B6/rij6);
                Ep += Ep_LJ;
            }
        }
    }

    cout << "Energi potensial dari molekul metana ruah: " \
         << Ep/(float)N << "kJ/mol";

    // hasil iterasi dimasukkan ke dalam file xyz
    ofstream file;
    file.open("metana-hitung_ep.xyz");
    file << N << "\n" << endl;

    int Nw = N/5;

    // iterasi dalam output array
    for(int m = 0; m < Nw; m++){
        file << setw(5) << "C" << setw(5) << " " \
             << fixed << setprecision(5) << Cx[m] << setw(5) << " " \
             << fixed << setprecision(5) << Cy[m] << setw(5) << " " \
             << fixed << setprecision(5) << Cz[m] << "\n";

        file << setw(5) << "H" << setw(5) << " " \
             << fixed << setprecision(5) << H1x[m] << setw(5) << " " \
             << fixed << setprecision(5) << H1y[m] << setw(5) << " " \
             << fixed << setprecision(5) << H1z[m] << "\n";

        file << setw(5) << "H" << setw(5) << " " \
             << fixed << setprecision(5) << H2x[m] << setw(5) << " " \
             << fixed << setprecision(5) << H2y[m] << setw(5) << " " \
             << fixed << setprecision(5) << H2z[m] << "\n";

        file << setw(5) << "H" << setw(5) << " " \
             << fixed << setprecision(5) << H3x[m] << setw(5) << " " \
             << fixed << setprecision(5) << H3y[m] << setw(5) << " " \
             << fixed << setprecision(5) << H3z[m] << "\n";

        file << setw(5) << "H" << setw(5) << " " \
             << fixed << setprecision(5) << H4x[m] << setw(5) << " " \
             << fixed << setprecision(5) << H4y[m] << setw(5) << " " \
             << fixed << setprecision(5) << H4z[m] << "\n";
    }   

    file.close();

    return 0;
}
