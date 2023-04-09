#include <iostream>
#include <vector>
#include <Eigen/Dense> // inclure la librairie Eigen
#include <fstream>  //permet de g�rer des flux d'entr�/sortie avec des fichiers externes

using namespace Eigen;
using namespace std;

int transient() {
    cout << "In transient mode." << endl;
    double l_x{ 1 };  //longueur de la barre en [m]
    int n_x{ 2 }; //nb de volumes de contr�le selon x
    double T_W{ 20 };  //conditions aux limites en W en [C�]
    double T_E{ 20 };  //conditions aux limites en E en [C�]
    double R{ 0.05 };  //diffusivit� thermique R en [m^2/s] (R = Lambda/(rho*c)
    double t_final{ 1 }; //dur�e de la simulation en [s]
    int n_t{ 10 };  //nb de pas de temps
    double T_init{ 30 };  //temp�rature initiale de la barre

    //flux d'entr� depuis un fichier
    string const fichier_inputs("E:\\Visual_Studio_Projects\\Physique\\solver_thermique_test\\x64\\Debug\\inputs_transient.txt");
    ifstream mon_flux_entre(fichier_inputs.c_str());
    string ligne;
    string nom_variable;
    double valeur_variable;
    if (mon_flux_entre) {
        while (getline(mon_flux_entre, ligne))
        {
            mon_flux_entre >> nom_variable;
            mon_flux_entre >> valeur_variable;
            if (nom_variable == "l_x") {
                l_x = valeur_variable;
            }
            else if (nom_variable == "n_x") {
                n_x = valeur_variable;
            }
            else if (nom_variable == "T_W") {
                T_W = valeur_variable;
            }
            else if (nom_variable == "T_E") {
                T_E = valeur_variable;
            }
            else if (nom_variable == "R") {
                R = valeur_variable;
            }
            else if (nom_variable == "t_final") {
                t_final = valeur_variable;
            }
            else if (nom_variable == "n_t") {
                n_t = valeur_variable;
            }
            else if (nom_variable == "T_0") {
                T_init = valeur_variable;
            }
        }
    }
    else {
        cout << "ERREUR: Impossible d'ouvrir le fichier '" << fichier_inputs << "' en lecture." << endl;
    }
    mon_flux_entre.close();

    //def mesh spatial
    double dx{ l_x / static_cast<double>(n_x) };
    /*vector<double> x(n_x);  //x est seulement pour faire de la visualisation
    for (int i = 0; i < x.size(); i++) {
        x[i] = dx/2 + i*dx;
        cout << "x[i] = " << x[i] << endl;
    }*/

    //def mesh temporel
    double dt{ t_final / static_cast<double>(n_t) };
    VectorXd t(n_t);
    t.setZero();
    for (int i = 0; i < n_t; i++) {
        t(i) = (i+1) * dt;
    }
    cout << "t = " << endl << t << endl;

    VectorXd T_0(n_x);
    T_0.setOnes();
    T_0 *= T_init;
    cout << "T_0 = " << endl << T_0 << endl;

    MatrixXd T(n_x, n_t + 1);
    T.setZero();
    for (int i = 0; i < n_x; i++) {
        T(i, 0) = T_0(i);
    }
    cout << "T = " << endl << T << endl;


    //materiaux infos
    double R_w{ R };
    double R_e{ R };
    double h{ 1 };  //l'�paisseur est �gale � l'unit�
    double dy{ 1 };  //la largeur est �gale � l'unit�
    double A_w{ dy * h };  //Section transversalle
    double A_e{ dy * h };
    double dx_w{ dx };  //Distance entre les points X et W
    double dx_e{ dx };  //Distance entre les points X et E
    double a_W{ R_w * A_w / dx_w };  //Coef devant les thermes d'int�r�ts
    double a_E{ R_e * A_e / dx_e };  //Coef devant les thermes d'int�r�ts
    double a_P0{ dx * dy * h / dt };
    double a_P{ a_P0 + a_E + a_W };
    cout << "l_x = " << l_x << endl;
    cout << "n_x = " << n_x << endl;
    cout << "T_W = " << T_W << endl;
    cout << "T_E = " << T_E << endl;
    cout << "R = " << R << endl;
    cout << "a_W = " << a_W << endl;
    cout << "a_E = " << a_E << endl;
    cout << "a_P0 = " << a_P0 << endl;
    cout << "a_P = " << a_P << endl;

    /*Ecriture des 3 types d'�quations
    -----------------
    | G |   C   | D |
    -----------------
    eqa(G)  = -0 + a_P * T(i, t+1) - aE * T(i+1, t+1) = aW * TW + a_P0 * T(i, t)
    eqa(C)  = -aW * T(i-1, t+1) + a_P * T(i, t+1) - aE * T(i+1, t+1) = a_P0 * T(i, t)
    eqa(D)  = -aW * T(i-1, t+1) + a_P * T(i, t+1) - 0 = a_P0 * T(i, t) + aE * TE
    */

    //formation de la matrice[A]
    MatrixXd A(n_x, n_x);
    A.setZero();
    A(0, 0) = a_P;
    A(0, 1) = -a_E;
    for (int i = 1; i < n_x - 1; i++) {
        A(i, i - 1) = -a_W;
        A(i, i) = a_P;
        A(i, i + 1) = -a_E;
    }
    A(n_x - 1, n_x - 2) = -a_W;
    A(n_x - 1, n_x - 1) = a_P;
    cout << "A = " << endl << A << endl;

    FullPivLU<MatrixXd> lu(A); // D�composition LU de A

    //formation du vecteur [b]
    MatrixXd b(n_x, n_t);
    b.setOnes();
    for (int t = 0; t < n_t; t++) {
        b(0, t) = a_W * T_W + a_P0 * T(0, t);
        for (int x = 1; x < n_x-1; x++) {
            b(x, t) = a_P0 * T(x, t);
        }
        b(n_x-1, t) = a_P0 * T(n_x-1, t) + a_E * T_E;
        //cout << "b = " << endl << b.col(t) << endl;

        //r�solution du syst�me d'�quation [A][T]=[b]
        T.col(t + 1) = lu.solve(b.col(t));
    }

    cout << "b = " << endl << b << endl;
    cout << "T = " << endl << T << endl;

    VectorXd T_W_vec(n_t + 1);
    T_W_vec.setOnes();
    T_W_vec *= T_W;

    VectorXd T_E_vec(n_t + 1);
    T_E_vec.setOnes();
    T_E_vec *= T_E;

    MatrixXd temperatures(n_x + 2, n_t + 1);
    temperatures.setZero();
    temperatures.row(0) = T_W_vec;
    for (int i = 1; i < n_x + 1; i++) {
        for (int j = 0; j < n_t+1; j++) {
            temperatures(i, j) = T(i-1, j);
        }
    }
    temperatures.row(n_x + 1) = T_E_vec;

    cout << "Solution du systeme lineaire A*t=b :" << endl;
    cout << "temperatures = " << endl << temperatures << endl;

    return 0;
}