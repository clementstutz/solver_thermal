#include "tools.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense> // inclure la librairie Eigen
#include <fstream>  //permet de gérer des flux d'entré/sortie avec des fichiers externes
#include <chrono>


using namespace Eigen;
using namespace std;

int reading_input(string const& fichier_input, double& l_x, int& n_x, double& T_W, double& T_E, double& R) {
    //flux d'entré depuis un fichier
    ifstream flux_entre(fichier_input.c_str());
    string ligne;
    string nom_variable;
    double valeur_variable;
    if (!flux_entre) {
        cerr << "ERREUR: Impossible d'ouvrir le fichier '" << fichier_input << "' en lecture." << endl;
        return 1;
    }

    while (getline(flux_entre, ligne))
    {
        flux_entre >> nom_variable;
        flux_entre >> valeur_variable;
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
    
    flux_entre.close();
    return 0;
}

int writing_output(string const& fichier_output, double& l_x, int& n_x, double& T_W, double& T_E, double& R, VectorXd& x, VectorXd& T) {
    // recuperation de la date
    string current_time = get_current_time_formatted();

    //flux de sorti vers un fichier
    //string const fichier_output("output.out");
    ofstream flux_sortie(fichier_output.c_str());
    if (!flux_sortie) {
        cerr << "ERREUR: Impossible d'ouvrir le fichier '" << fichier_output << "' en ecriture." << endl;
        return 1;
    }

    flux_sortie << "Date : " << current_time << endl;
    flux_sortie << "" << endl;
    flux_sortie << "Input data used for the calculation :" << endl;
    flux_sortie << "l_x " << l_x << endl;
    flux_sortie << "n_x " << n_x << endl;
    flux_sortie << "T_W " << T_W << endl;
    flux_sortie << "T_E " << T_E << endl;
    flux_sortie << "R " << R << endl;
    flux_sortie << "" << endl;
    flux_sortie << "Output of the calculation :" << endl;
    flux_sortie << "x =" << endl << x << endl;
    flux_sortie << "T = " << endl << T << endl;

    flux_sortie.close();

    return 0;
}

int stationary(string const& fichier_input, string const& fichier_output) {
    cout << "In stationary mode." << endl;
    double l_x{ 1 };  //longueur de la barre en [m]
    int n_x{ 2 }; //nb de volumes de contrôle selon x
    double T_W{ 20 };  //conditions aux limites en W en [C°]
    double T_E{ 20 };  //conditions aux limites en E en [C°]
    double R{ 0.05 };  //diffusivité thermique R en [m^2/s] (R = Lambda/(rho*c)

    //reading the input data
    int RETURN_reading_input = reading_input(fichier_input, l_x, n_x, T_W, T_E, R);
    if (RETURN_reading_input != 0) {
        return 1;
    }

    //def spatial mesh
    double dx{ l_x / static_cast<double>(n_x) };
    VectorXd x(n_x + 2);  //x est seulement pour faire de la visualisation
    x.setZero();
    for (int i = 1; i < n_x + 1; i++) {
        x(i) = dx / 2 + (i - 1) * dx;
    }
    x(n_x + 1) = l_x;
    //cout << "x = " << endl << x << endl;

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
    //cout << "A = " << endl << A << endl;

    //formation du vecteur [b]
    VectorXd b(n_x);
    b.setZero();
    b(0) = -2 * a_W * T_W;
    b(n_x - 1) = -2 * a_E * T_E;
    //cout << "b = " << endl << b << endl;

    //résolution du système d'équation [A][T]=[b]
    // Inversion de A
    MatrixXd A_inv(n_x, n_x);
    A_inv = A.inverse();
    VectorXd T(n_x);
    //auto start_time = chrono::high_resolution_clock::now();
    T = A_inv * b;
    //auto end_time = chrono::high_resolution_clock::now();
    //auto elapsed_time_ms = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();

    cout << "Solution du systeme lineaire A*t=b :" << endl;
    cout << "T = " << endl << T << endl;

    int RETURN_writing_output = writing_output(fichier_output, l_x, n_x, T_W, T_E, R, x, T);


    return 0;
}
