#include "tools.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense> // inclure la librairie Eigen
#include <fstream>  //permet de gérer des flux d'entré/sortie avec des fichiers externes
#include <chrono>


using namespace Eigen;
using namespace std;

int reading_input(string const& fichier_input, double& l_x, int& n_x, double& T_W, double& T_E, double& R, double& t_final, int& n_t, double& T_init, double& taux_convergence) {
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
        else if (nom_variable == "t_final") {
            t_final = valeur_variable;
        }
        else if (nom_variable == "n_t") {
            n_t = valeur_variable;
        }
        else if (nom_variable == "T_0") {
            T_init = valeur_variable;
        }
        else if (nom_variable == "taux_convergence") {
            taux_convergence = valeur_variable;
        }
    }
    
    flux_entre.close();
    return 0;
}

int writing_output(string const& fichier_output, double& l_x, int& n_x, double& T_W, double& T_E, double& R, double& t_final, int& n_t, double& T_init, VectorXd& x, MatrixXd& T, VectorXd& t) {
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
    flux_sortie << "t_final " << t_final << endl;
    flux_sortie << "n_t " << n_t << endl;
    flux_sortie << "T_init " << T_init << endl;
    flux_sortie << "" << endl;
    flux_sortie << "Output of the calculation :" << endl;
    flux_sortie << "x =" << endl << x << endl;
    flux_sortie << "t =" << endl << t << endl;
    flux_sortie << "T = " << endl << T << endl;

    flux_sortie.close();

    return 0;
}

double error_global(VectorXd& v1, VectorXd& v2) {

    double error{ 0 };

    for (int i = 0; i < v1.size(); i++) {
        error += (v2(i) - v1(i)) / v1(i);
    }

    return error /= v1.size();
}

int add_vec_to_mat(vector<vector<double>>& matrice, int& new_vec_size) {
    // Ajoute un vecteur supplementaire au vecteur de vecteur matrice
    matrice.resize(matrice.size() + 1);

    //definition de la taille du vecteur ajoute
    matrice[matrice.size() - 1].resize(new_vec_size);

    //remplissage du nouveau vecteur
    //for (int i = 0; i < new_vec_size; i++) {
    //    //matrice_T[matrice_T_indice_size][i] = values[i];

    //    //afichage du contenu du nouveau vecteur
    //    //cout << matrice[matrice.size() - 1][i] << endl;
    //}

    return 0;
}

int transient(string const& fichier_input, string const& fichier_output) {
    cout << "In transient mode." << endl;
    double l_x{ 1 };  //longueur de la barre en [m]
    int n_x{ 2 }; //nb de volumes de contrôle selon x
    double T_W{ 20 };  //conditions aux limites en W en [C°]
    double T_E{ 20 };  //conditions aux limites en E en [C°]
    double R{ 0.05 };  //diffusivité thermique R en [m^2/s] (R = Lambda/(rho*c)
    double t_final{ 1 }; //durée de la simulation en [s]
    int n_t{ 10 };  //nb de pas de temps
    double T_init{ 30 };  //température initiale de la barre
    double taux_convergence{ 0.05 };

    //reading the input data
    int RETURN_reading_input = reading_input(fichier_input, l_x, n_x, T_W, T_E, R, t_final, n_t, T_init, taux_convergence);
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

    //def temporal mesh
    double dt{ t_final / static_cast<double>(n_t) };
    VectorXd t(n_t+1);
    t.setZero();
    for (int i = 0; i < n_t+1; i++) {
        t(i) = (i) * dt;
    }
    //cout << "t = " << endl << t << endl;

    VectorXd T_0(n_x);
    T_0.setOnes();
    T_0 *= T_init;
    //cout << "T_0 = " << endl << T_0 << endl;

    MatrixXd T(n_x, n_t + 1);
    T.setZero();
    for (int i = 0; i < n_x; i++) {
        T(i, 0) = T_0(i);
    }
    cout << "T = " << endl << T << endl;

    vector<vector<double>> matrice_T(1); // crer un vecteur de vecteurs
    matrice_T[matrice_T.size() - 1].resize(n_x);
    cout << "matrice_T[" << matrice_T.size() - 1 << "][i] = " << endl;
    for (int i = 0; i < n_x; i++) {
        matrice_T[matrice_T.size() - 1][i] = T_0(i);
        cout << matrice_T[matrice_T.size() - 1][i] << endl;
    }

    

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

    /*Ecriture des 3 types d'équations
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
    //cout << "A = " << endl << A << endl;

    // Décomposition LU de A
    FullPivLU<MatrixXd> lu(A);

    //formation du vecteur [b]
    MatrixXd b(n_x, n_t);
    b.setOnes();
    for (int t = 0; t < n_t; t++) {
        b(0, t) = a_W * T_W + a_P0 * T(0, t);
        for (int x = 1; x < n_x - 1; x++) {
            b(x, t) = a_P0 * T(x, t);
        }
        b(n_x - 1, t) = a_P0 * T(n_x - 1, t) + a_E * T_E;
        //cout << "b = " << endl << b.col(t) << endl;

        //résolution du système d'équation [A][T]=[b]
        add_vec_to_mat(matrice_T, n_x);
        
        VectorXd T1(n_x);
        T1 = lu.solve(b.col(t));
        
        cout << "matrice_T[" << matrice_T.size() - 1 << "][i] = " << endl;
        for (int i = 0; i < n_x; i++) {
            matrice_T[matrice_T.size() - 1][i] = T1(i);
            cout << matrice_T[matrice_T.size() - 1][i] << endl;
        }


        T.col(t + 1) = lu.solve(b.col(t));

        VectorXd T_1(n_x);
        T_1 = T.col(t);

        VectorXd T_2(n_x);
        T_2 = T.col(t + 1);

        
        double error{ 0 };
        error = error_global(T_1, T_2);
        if (error <= taux_convergence) {
            cout << "has converged. error = " << error << endl;
            break;
        }
        cout << "has NOT converged. error = " << error << endl;
    }
    //cout << "b = " << endl << b << endl;
    cout << "T = " << endl << T << endl;

    VectorXd T_W_vec(n_t + 1);
    T_W_vec.setOnes();
    T_W_vec *= T_W;

    VectorXd T_E_vec(n_t + 1);
    T_E_vec.setOnes();
    T_E_vec *= T_E;

    MatrixXd TT(n_x + 2, n_t + 1);
    TT.setZero();
    TT.row(0) = T_W_vec;
    for (int i = 1; i < n_x + 1; i++) {
        for (int j = 0; j < n_t + 1; j++) {
            TT(i, j) = T(i - 1, j);
        }
    }
    TT.row(n_x + 1) = T_E_vec;

    cout << "Solution du systeme lineaire A*t=b :" << endl;
    cout << "T = " << endl << TT << endl;

    int RETURN_writing_output = writing_output(fichier_output, l_x, n_x, T_W, T_E, R, t_final,  n_t,  T_init, x, T, t);

    return 0;
}