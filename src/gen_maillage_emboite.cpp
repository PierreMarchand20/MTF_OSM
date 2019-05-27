#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;
const double pi = 3.141592653589793;

int main(int argc, char* argv[]){

    // Check the number of parameters
    if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " nombre d'interface " << std::endl;
        /* "Usage messages" are a conventional way of telling the user
        * how to run a program if they enter the command incorrectly.
        */
        return 1;
    }
    int nc = std::stoi(argv[1]);


    double h = 0.1;
    double ratio = 1.;
    int nb_node = 0;
    vector<double> rad(nc),theta(nc);
    vector<int> nv(nc);

    for(int j=0; j<nc; j++){
        rad[j]   = ratio*(nc-j);
        nv[j]    = floor(rad[j]*2*pi/h);
        theta[j] = 2*pi/double(nv[j]);
        nb_node += nv[j];
    }

    ofstream meshfile;
    meshfile.open(("../meshes/maillage_emboite_"+std::to_string(nc)+".msh").c_str());
    meshfile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
    meshfile << nb_node << endl;
    int num_node = 1;
    for(int j = 0; j<nc; j++){
        double r,x,y,z,t,dt;
        r  = rad[j];
        dt = theta[j];

        for(int k=0; k<nv[j]; k++){
            t = k*dt;
            x = r*cos(t);
            y = r*sin(t);
            z = 0.;

            meshfile << left << setw(10) << num_node;
            meshfile << left << setw(15) << x;
            meshfile << left << setw(15) << y;
            meshfile << left << setw(15) << z << endl;
            num_node++;
        }
    }

    meshfile << "$EndNodes\n";
    meshfile << "$Elements\n";
    meshfile << nb_node << endl;
    int num_elt = 1;
    for(int j = 0; j<nc; j++){
        int nbv = nv[j];
        for(int k=0; k<nv[j]; k++){
            meshfile << left << setw(6) << num_elt+k << "1 2 ";
            meshfile << j << " 1"<< "\t";
            meshfile << num_elt+ k%nbv     << "\t";
            meshfile << num_elt+ (k+1)%nbv << "\n";
        }
        num_elt += nv[j];
    }
    meshfile << "$EndElements" << endl;
    meshfile.close();

}
