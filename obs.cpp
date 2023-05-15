#include <iostream>
#include <cmath>
using namespace std;

int writefile(int a, int b, int r);
float dist(int ii, int jj);

#define CENTERX 64
#define CENTERY 64

int writefile(int a, int b, int r){
    string filename = "obstacles_" + to_string(a) + "x" + to_string(b) + ".dat";
    FILE* fp = fopen(filename.c_str(), "w");
        
    for (int jj = 0 ; jj < b; jj += b-1){
        for (int ii = 0; ii < a; ii ++){
            fprintf(fp, "%d %d %d\n", ii, jj, 1);
        }
    }

    for (int ii = CENTERX - r; ii <= CENTERX+r; ii++){
        for (int jj = CENTERY - r; jj <= CENTERY + r; jj ++){
            if (dist(ii, jj) < r){
                fprintf(fp, "%d %d %d\n", ii, jj, 1);
            }
        }
    }
    fclose(fp);
    return EXIT_SUCCESS;
}

float dist(int ii, int jj){
    return sqrt(pow(ii-CENTERX, 2) + pow(jj-CENTERY, 2));
}

int main(){
    int a = 512;
    int b = 256;
    int r = 3;
    writefile(a, b, r);
}
