#include <iostream>
#include <chrono>

using namespace std;

string get_current_time_formatted() {
    time_t t = time(NULL);
    struct tm tm;

    if (localtime_s(&tm, &t) != 0) {
        cerr << "Erreur lors de la récupération de la date et de l'heure.\n";
        return "";
    }

    stringstream ss;
    ss << put_time(&tm, "%d/%m/%Y : %H:%M:%S");

    return ss.str();
}