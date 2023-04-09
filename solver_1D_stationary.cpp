#include <iostream>
#include <vector>
#include <Eigen/Dense> // inclure la librairie Eigen
#include <fstream>  //permet de gérer des flux d'entré/sortie avec des fichiers externes

using namespace Eigen;
using namespace std;

int stationary() {
    cout << "In stationary mode." << endl;
    double l_x{ 1 };  //longueur de la barre en [m]
    int n_x{ 2 }; //nb de volumes de contrôle selon x
    double T_W{ 20 };  //conditions aux limites en W en [C°]
    double T_E{ 20 };  //conditions aux limites en E en [C°]
    double R{ 0.05 };  //diffusivité thermique R en [m^2/s] (R = Lambda/(rho*c)

    //flux d'entré depuis un fichier
    string const fichier_inputs("E:\\Visual_Studio_Projects\\Physique\\solver_thermique_test\\x64\\Debug\\inputs_stationary.txt");
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
        }
    }
    else {
        cout << "ERREUR: Impossible d'ouvrir le fichier '" << fichier_inputs << "' en lecture." << endl;
    }
    mon_flux_entre.close();

    //def mesh
    double dx{ l_x / static_cast<double>(n_x) };
    /*vector<double> x(n_x);  //x est seulement pour faire de la visualisation
    for (int i = 0; i < x.size(); i++) {
        x[i] = dx/2 + i*dx;
        cout << "x[i] = " << x[i] << endl;
    }*/

    //materiaux infos
    double R_w{ R };
    double R_e{ R };
    double h{ 1 };  //l'épaisseur est égale à l'unité
    double dy{ 1 };  //la largeur est égale à l'unité
    double A_w{ dy * h };  //Section transversalle
    double A_e{ dy * h };
    double dx_w{ dx };  //Distance entre les points X et W
    double dx_e{ dx };  //Distance entre les points X et E
    double a_W{ R_w * A_w / dx_w };  //Coef devant les thermes d'intérêts
    double a_E{ R_e * A_e / dx_e };  //Coef devant les thermes d'intérêts
    cout << "l_x = " << l_x << endl;
    cout << "n_x = " << n_x << endl;
    cout << "T_W = " << T_W << endl;
    cout << "T_E = " << T_E << endl;
    cout << "R = " << R << endl;
    cout << "a_W = " << a_W << endl;
    cout << "a_E = " << a_E << endl;

    /*Ecriture des 3 types d'équations
    -----------------
    | G |   C   | D |
    -----------------
    eqa(G) = 0 - (2 * a_W + a_E) * T(1) + a_E * T(2) = -2 * a_W * T_W
    eqa(C) = a_W * T(i - 1) - (a_W + a_E) * T(i) + a_E * T(i + 1) = 0
    eqa(D) = a_W * T(n_x - 1) - (a_W + 2 * a_E) * T(n_x) + 0 = -2 * a_E * T_E
    */

    //formation de la matrice[A]
    MatrixXd A(n_x, n_x);
    A.setZero();
    A(0, 0) = -(2 * a_W + a_E);
    A(0, 1) = a_E;
    for (int i = 1; i < n_x - 1; i++) {
        A(i, i - 1) = a_W;
        A(i, i) = -(a_W + a_E);
        A(i, i + 1) = a_E;
    }
    A(n_x - 1, n_x - 2) = a_W;
    A(n_x - 1, n_x - 1) = -(a_W + 2 * a_E);
    cout << "A = " << endl << A << endl;

    //formation du vecteur [b]
    VectorXd b(n_x);
    b.setZero();
    b(0) = -2 * a_W * T_W;
    b(n_x - 1) = -2 * a_E * T_E;
    cout << "b = " << endl << b << endl;

    //résolution du système d'équation [A][T]=[b]
    VectorXd t(n_x);
    t = A.fullPivLu().solve(b); // résoudre le système linéaire

    cout << "Solution du systeme lineaire A*t=b :" << endl;
    cout << "T = " << endl << t << endl; // afficher la solution

    return 0;
}
