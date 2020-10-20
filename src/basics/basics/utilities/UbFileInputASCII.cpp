#include <basics/utilities/UbFileInputASCII.h>
#include <cstring>

using namespace std;

UbFileInputASCII::UbFileInputASCII(string filename)
{
    this->filename         = filename;
    this->commentindicator = 'C';

    infile.open(filename.c_str());

    // if(!infile) UB_THROW( UbException((string)("UbFileInputASCII::UbFileInputASCII(string filename, int how) couldn't
    // open file:\n "+filename)) );
}
/*==========================================================*/
bool UbFileInputASCII::open(string filename)
{
    infile.close();
    infile.clear(); // setzt flags zurueck

    this->filename = filename;
    infile.open(this->filename.c_str());

    return infile.is_open();
}
/*==========================================================*/
int UbFileInputASCII::readInteger()
{
    int dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
long long UbFileInputASCII::readLongLong()
{
    long long dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
string UbFileInputASCII::getFileName() { return this->filename; }

/*==========================================================*/
void UbFileInputASCII::skipLine()
{
    string dummy;
    getline(infile, dummy);
}
/*==========================================================*/
void UbFileInputASCII::readLine()
{
    string dummy;
    getline(infile, dummy);
}
/*==========================================================*/
string UbFileInputASCII::readStringLine()
{
    string dummy;
    getline(infile, dummy);
    return dummy;
}
/*==========================================================*/
string UbFileInputASCII::readLineTill(char stop)
{
    string dummy;
    getline(infile, dummy, stop);
    return dummy;
}
/*==========================================================*/
string UbFileInputASCII::parseString()
{
    string dummy;
    getline(infile, dummy, ' ');
    return dummy;
}
/*==========================================================*/
double UbFileInputASCII::readDouble()
{
    double dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
float UbFileInputASCII::readFloat()
{
    float dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
string UbFileInputASCII::readString()
{
    string dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
char UbFileInputASCII::readChar()
{
    int dummy;
    infile >> dummy;
    return (char)dummy;
}
/*==========================================================*/
std::size_t UbFileInputASCII::readSize_t()
{
    std::size_t dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
void UbFileInputASCII::setPosAfterLineWithString(const string &var)
{
    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen
    char line[512];
    do {
        infile.getline(line, 512);
        if (infile.eof())
            UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> string " + var +
                                                " wasn't found in " + this->filename));
    } while (strstr(line, var.c_str()) != line); // Ende Schleife, wenn varname ganz in zeile vorkommt
}
/*==========================================================*/
bool UbFileInputASCII::containsString(const string &var)
{
    infile.clear(); // setzt den EOF-Status zurueck (wird durch infile.seekg() NICHT getan!!!)

    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen
    char line[512];
    do {
        infile.getline(line, 512);
        if (infile.eof())
            return false;
    } while (strstr(line, var.c_str()) != line); // Ende Schleife, wenn varname ganz in zeile vorkommt

    return true;
}
/*==========================================================*/
int UbFileInputASCII::readIntegerAfterString(const string &var)
// last change [10.3.2004] at [9:46]
// suchts in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
// z.B. timesteps 9
{
    infile.clear(); // setzt den EOF-Status zurueck (wird durch infile.seekg() NICHT getan!!!)

    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen

    char line[512];

    do {
        infile.getline(line, 512);
        if (infile.eof())
            UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> " + var +
                                                " wasn't found in " + this->filename));
    } while (strstr(line, var.c_str()) != line); // Ende Schleife, wenn varname ganz in zeile vorkommt

    strcpy(line, (line + strlen(var.c_str()))); // zeile um "varname" kuerzen
    while ((line[0] == ' ') || (line[0] == '\t'))
        strcpy(line, (line + 1)); // Whitespaces entfernen

    return (atoi(line)); // Umwandlung in int
}
/*==========================================================*/
// last change [10.3.2004] at [9:46]
// sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
// z.B. nue 9.5
double UbFileInputASCII::readDoubleAfterString(const string &var)
{
    infile.clear(); // setzt den EOF-Status zurueck (wird durch infile.seekg() NICHT getan!!!)

    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen

    char line[512];

    do {
        infile.getline(line, 512);
        if (infile.eof())
            UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> " + var +
                                                " wasn't found in " + this->filename));
    } while (/*!strncmp(varname,line,sizeof(varname))==0*/ strstr(line, var.c_str()) !=
             line); // Ende Schleife, wenn varname ganz in zeile vorkommt

    strcpy(line, (line + strlen(var.c_str()))); // zeile um "varname" kuerzen
    while ((line[0] == ' ') || (line[0] == '\t'))
        strcpy(line, (line + 1)); // Whitespaces entfernen

    return (atof(line)); // Umwandlung in double
}
/*==========================================================*/
//  [9.9.2002]
// liefert string-Wert der hinter dem uebergebenen char feld in der datei infile steht
// zudem wird der wert in die uebergebene variable value uebertragen (falls man das ergebniss als char benoetig)
string UbFileInputASCII::readStringAfterString(const string &var) //,char *value)
{
    infile.clear(); // setzt den EOF-Status zurueck (wird durch infile.seekg() NICHT getan!!!)

    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen

    char line[512];
    // string line_copy[512];

    do {
        infile.getline(line, 512);
        if (infile.eof())
            UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> " + var +
                                                " wasn't found in " + this->filename));
    } while (strstr(line, var.c_str()) != line); // Ende Schleife, wenn varname ganz in zeile vorkommt

    strcpy(line, (line + strlen(var.c_str()))); // zeile um "varname" kuerzen
    while ((line[0] == ' ') || (line[0] == '\t'))
        strcpy(line, (line + 1)); // Whitespaces entfernen

    char *p;
    p = strtok(line, " ");  // schneidet alles "ab und inklusive space " nach namen ab
    p = strtok(line, "\t"); // schneidet alles "ab und inklusive tab   " nach namen ab

    return static_cast<string>(p); // Umwandlung in string
}
/*==========================================================*/
// last change [10.3.2004] at [9:46]
// sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
// z.B. nue 9.5
bool UbFileInputASCII::readBoolAfterString(const string &var)
{
    if (this->readStringAfterString(var) == "true")
        return true;
    else if (this->readStringAfterString(var) == "false")
        return false;
    else
        UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> expression after " + var +
                                            " is not equal to 'true' or 'false' in " + this->filename));
}
/*==========================================================*/
// last change [10.3.2004] at [9:46]
// sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
// z.B. nue 9.5
bool UbFileInputASCII::readBool()
{
    string tmp = this->readString();
    if (tmp == "true")
        return true;
    else if (tmp == "false")
        return false;
    else
        UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> expression=\"" + tmp +
                                            "\" is not equal to 'true' or 'false' in " + this->filename));
}
