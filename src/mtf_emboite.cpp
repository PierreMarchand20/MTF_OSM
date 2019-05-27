#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/mat_struct.hpp>

using namespace bemtool;
using namespace std;


Cplx inc_wave(const Real& k, const R3& x){
    R3 d; d[0]=1.; return exp( iu*k*(d,x) );
}

C3 grad_inc_wave(const Real& k, const R3& x){
    R3 d; d[0]=1.; return iu*k*exp( iu*k*(d,x) )*d;
}

int main(int argc, char* argv[]){

    // Check the number of parameters
    if (argc < 3) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " nombre d'interface \t type" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
        * how to run a program if they enter the command incorrectly.
        */
        return 1;
    }

    // Data
    string nc_str= argv[1];
    int nc = std::stoi(argv[1]);
    string type_str = argv[2];
    int type = std::stoi(argv[2]);
    std::vector<double> kappa(nc+1),mu(nc+1);
    srand (1);
    for (int i=0;i<kappa.size();i++){

        if (type==0){
            kappa[i] =1;
            double pert = 2*((double) rand() / (double)(RAND_MAX));
            kappa[i]+=pert;
            mu[i] = 1;
        }
        else if (type==1){
            kappa[i] = 1;
            mu[i] = 1;
        }

        else if (type==2){
            kappa[i] = (i%2==0) ? 1 : 3;
            mu[i] = 1;
        }
        else if (type==3){
            kappa[i] =1;
            double pert = 40*((double) rand() / (double)(RAND_MAX));
            kappa[i]+=pert;
            mu[i] = 1;
        }
        else if (type==4){
            kappa[i] = 1;
            mu[i] = 1;
            double pert = 100*((double) rand() / (double)(RAND_MAX));
            mu[i]+=pert;
        }     
        else if (type==5){
            kappa[i] = (i%2==0) ? 1 : 10;
            mu[i] = 1;
        }
        std::cout << kappa[i] << " "<<mu[i]<<endl;
    }

    // Noeuds du maillage
    Geometry node(("../meshes/maillage_emboite_"+nc_str+".msh").c_str());
    cout << "nb node: " << NbNode(node) << endl;

    const int NbInt = nc;
    const int NbDom = NbInt+1;

    std::vector<Mesh1D> Gamma(NbInt);
    for(int j=0; j<NbInt; j++){
        Gamma[j].Load(node,j);
    }

    std::vector<Mesh1D> Omega(NbDom);
    Omega[0]  = unbounded;
    Omega[0] += Gamma[0];
    for(int j=0; j<NbDom-2; j++){
        Omega[j+1] += Gamma[j];
        Omega[j+1] += Gamma[j+1];
    }
    Omega[NbDom-1] += Gamma[NbInt-1];

    for(int j=0; j<NbDom; j++){
        Orienting(Omega[j]);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    Degres de liberte     //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    typedef Dof<P1_1D> DofType;
    std::vector<DofType> dof;
    for (int i =0;i<NbDom;i++){
        dof.emplace_back(Omega[i]);
    }
    // DofType dof[NbDom] = { DofType(Omega[0]),
    //     DofType(Omega[1]),
    //     DofType(Omega[2]),
    //     DofType(Omega[3]),
    //     DofType(Omega[4]),
    //     DofType(Omega[5]) };

    int nb_dof[NbDom];
    for(int j=0; j<NbDom; j++){
        nb_dof[j]=NbDof(dof[j]);
        cout << "nb_dof[" << j << "] = " << nb_dof[j] << endl;
    }

    int nb_dof_tot = 0;
    int nD[NbDom], nN[NbDom];
    for(int j=0; j<NbDom; j++){
        nD[j] = nb_dof_tot;
        nN[j] = nD[j]+nb_dof[j];
        nb_dof_tot += 2*nb_dof[j];
    }

    cout << "nb_dof_tot = " << nb_dof_tot << endl;


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    Assembly of matrices      //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DenseMatrix<Cplx>  A(nb_dof_tot,nb_dof_tot);
    DenseMatrix<Cplx>  M(nb_dof_tot,nb_dof_tot);
    SparseMatrix<Cplx>  Mass(nb_dof_tot,nb_dof_tot);


    for(int I=0; I<NbDom; I++){
        const Mesh1D&  mesh = Omega[I];
        const DofType& dofI = dof[I];
        const int ne        = NbElt(mesh);

        BIOp<HE_SL_2D_P1xP1 >  V(mesh,mesh,kappa[I]);
        BIOp<HE_DL_2D_P1xP1 >  K(mesh,mesh,kappa[I]);
        BIOp<HE_TDL_2D_P1xP1> TK(mesh,mesh,kappa[I]);
        BIOp<HE_HS_2D_P1xP1 >  W(mesh,mesh,kappa[I]);

        for(int j=0; j<ne; j++){
            for(int k=0; k<ne; k++){

                N2 jD = nD[I]+dofI[j]; 	N2 kD = nD[I]+dofI[k];
                N2 jN = nN[I]+dofI[j]; 	N2 kN = nN[I]+dofI[k];

                // Assemblage operateurs integraux
                A(jD,kD) += 2*(1./mu[I])*W (j,k);
                A(jN,kD) += 2*K (j,k);
                A(jD,kN) += 2*TK(j,k);
                A(jN,kN) += 2*mu[I]*V (j,k);

            }

            // Assemblage matrice de masse
            const Elt1D& e = Omega[I][j];
            M(nN[I]+dofI[j], nD[I]+dofI[j]) += MassP1(e);
            M(nD[I]+dofI[j], nN[I]+dofI[j]) += MassP1(e);

            Mass.Insert(nD[I]+dofI[j],nD[I]+dofI[j],MassP1(e));
            Mass.Insert(nN[I]+dofI[j],nN[I]+dofI[j],MassP1(e));

            // A(nN[I]+dofI[j], nD[I]+dofI[j]) += MassP1(e);
            // A(nD[I]+dofI[j], nN[I]+dofI[j]) += MassP1(e);
        }
    }

    Mass.Export(("../output/matrices/Mass_emboite_"+nc_str+"_"+type_str+".txt").c_str());


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    Operateur de transmission   //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SparseMatrix<Cplx>  P(nb_dof_tot,nb_dof_tot);

    N2 interface[NbInt];
    for(int j=0; j<NbInt; j++){
        interface[j][0]=j; interface[j][1]=j+1;
    }

    for(int I=0; I<NbInt; I++){
        const N2& J = interface[I];

        for(int j=0; j<NbElt(Gamma[I]); j++){
            const Elt1D& e = Gamma[I][j];
            int j0 = node[e][ Omega[J[0]] ];
            int j1 = node[e][ Omega[J[1]] ];

            N2 jD = nD[J[0]]+dof[J[0]][j0];
            N2 kD = nD[J[1]]+dof[J[1]][j1];
            N2 jN = nN[J[0]]+dof[J[0]][j0];
            N2 kN = nN[J[1]]+dof[J[1]][j1];

            P.Insert(jN,kD,(+1.)*MassP1(e));
            P.Insert(jD,kN,(-1.)*MassP1(e));
            P.Insert(kN,jD,(+1.)*MassP1(e));
            P.Insert(kD,jN,(-1.)*MassP1(e));

            // A(jN,kD) += (+1.)*MassP1(e);
            // A(jD,kN) += (-1.)*MassP1(e);
            // A(kN,jD) += (+1.)*MassP1(e);
            // A(kD,jN) += (-1.)*MassP1(e);

        }
    }

    P.Export(("../output/matrices/Prec_emboite_"+nc_str+"_"+type_str+".txt").c_str());

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //   Assemblage second membre  //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<Cplx> b(nb_dof_tot,0.);
    vector<Cplx> f(nb_dof_tot,0.);

    {
        const std::vector<R3>& nor = NormalTo(Omega[0]);
        for(int j=0; j<NbElt(Omega[0]); j++){
            const Real& k0 = kappa[0];
            const Elt1D& e = Omega[0][j];
            const N2& jD   = nD[0]+dof[0][j];
            const N2& jN   = nN[0]+dof[0][j];

            b[jD[0]] += 0.5*inc_wave(k0,e[0]);
            b[jD[1]] += 0.5*inc_wave(k0,e[1]);

            b[jN[0]] += 0.5*(grad_inc_wave(k0,e[0]),nor[j]);
            b[jN[1]] += 0.5*(grad_inc_wave(k0,e[1]),nor[j]);
        }
    }
    mv_prod(f,M,b);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    Calcul solution exacte   //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<Cplx> vex(nb_dof_tot,0.);
    vector<Cplx> uex(nb_dof_tot,0.);

    for(int I=0; I<NbDom; I++){
        const std::vector<R3>& nor = NormalTo(Omega[I]);
        for(int j=0; j<NbElt(Omega[I]); j++){
            const Real& k0 = kappa[I];
            const Elt1D& e = Omega[I][j];
            const N2& jD   = nD[I]+dof[I][j];
            const N2& jN   = nN[I]+dof[I][j];

            vex[jD[0]] += 0.5*inc_wave(k0,e[0]);
            vex[jD[1]] += 0.5*inc_wave(k0,e[1]);

            vex[jN[0]] += 0.5*(grad_inc_wave(k0,e[0]),nor[j]);
            vex[jN[1]] += 0.5*(grad_inc_wave(k0,e[1]),nor[j]);
        }
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    Export des resultats     //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ofstream file;

    file.open("../output/matrices/A_emboite_"+nc_str+"_"+type_str+".txt");
    for(int j=0; j<nb_dof_tot; j++){
        for(int k=0; k<nb_dof_tot; k++){
            Cplx z = A(j,k);
            file << z.real() << " ";
            file << z.imag() << " ";
        }
    }
    file.close();

    file.open("../output/matrices/M_emboite_"+nc_str+"_"+type_str+".txt");
    for(int j=0; j<nb_dof_tot; j++){
        for(int k=0; k<nb_dof_tot; k++){
                Cplx z = M(j,k);
                file << z.real() << " ";
        }
    }
    file.close();

    file.open("../output/matrices/P_emboite_"+nc_str+"_"+type_str+".txt");
    for(int j=0; j<nb_dof_tot; j++){
        for(int k=0; k<nb_dof_tot; k++){
            Cplx z = P(j,k);
            file << z.real() << " ";
        }
    }
    file.close();

    file.open("../output/matrices/f_emboite_"+nc_str+"_"+type_str+".txt");
    for(int j=0; j<f.size(); j++){
        Cplx z = f[j];
        file << z.real() << " ";
        file << z.imag() << " ";
    }
    file.close();

    file.open("../output/matrices/uex_emboite_"+nc_str+"_"+type_str+".txt");
    for(int j=0; j<f.size(); j++){
        Cplx z = vex[j];
        file << z.real() << " ";
        file << z.imag() << " ";
    }
    file.close();

}
