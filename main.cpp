#include <iostream>
#include "solver_1D_stationary.h"
#include "solver_1D_transient.h"
#include <fstream>
#include <iomanip>
#include <sstream>


using namespace std;

int main(int argc, char *argv[])
{
    // Vérification du nombre d'arguments
    if (argc < 3) {
        cerr << "Erreur : aucun fichier d'inputs ou d'outputs n'a ete specifie.\n";
        return 1;
    }

    // Récupération du chemin du fichier
    string fichier_inputs = argv[1];

    // Ouverture du fichier
    ifstream infile(fichier_inputs);
    if (!infile) {
        cerr << "Erreur : impossible d'ouvrir le fichier " << fichier_inputs << "\n";
        return 1;
    }

    // Récupération du chemin du fichier
    string fichier_outputs = argv[2];

    // Ouverture du fichier
    /*ifstream outfile(fichier_outputs);
    if (outfile) {
        cerr << "Erreur : le fichier de sortir " << fichier_inputs << " existe deja, veuillez le supprimer puis relancer le calcul.\n";
        return 1;
    }*/

    int RETURN_stationary = stationary(fichier_inputs, fichier_outputs);
    if (RETURN_stationary != 0) {
        return 1;
    }

    string const fichier_transient_inputs("E:\\Visual_Studio_Projects\\Physique\\solver_thermique_test\\x64\\Debug\\inputs_transient.txt");
    string const fichier_transient_outputs("E:\\Visual_Studio_Projects\\Physique\\solver_thermique_test\\x64\\Debug\\outputs_transient.out");
    int RETURN_transient = transient(fichier_transient_inputs, fichier_transient_outputs);
    if (RETURN_transient != 0) {
        return 1;
    }

    return 0;
}