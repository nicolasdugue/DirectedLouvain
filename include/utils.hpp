#ifndef UTILS_HPP
#define UTILS_HPP

#include <getopt.h>
#include <iostream>

using namespace std;

string filename = "";
string filename_part = "";
bool weighted = false;
double precision = 0.000001;
// No display done
int display_level = -1;
bool verbose = false;
bool reproducibility = false;
bool renumbering = true;

void usage(char * prog_name) {
    cerr << "usage: " << prog_name << " -f [-r] [-w] [-n] [-p partition_file] [-q epsilon] [-l display_level] [-v] [-h]" << endl;
    cerr << "usage: " << prog_name << " --file [--reproducibility] [--weigthed] [--renumbered] [--partition partition_file] [--precision epsilon] [--level display_level] --verbose] [--help]" << endl << endl;
    cerr << "-f|--file input file: file containing the graph to decompose in communities." << endl;
    cerr << "\tcan be any file under the edgelist format (any extension except \".bin\") _or_ a binary (\".bin\") file generated by the -r option" << endl;
    cerr << "[-r|--reproducibility] for the sake of fast reproducibility, the graph is written in binary form." << endl;
    cerr << "[-w|--weighted] \tread the graph as a weighted one (weights are set to 1 otherwise)." << endl;
    cerr << "[-n|--renumbered] \tto indicate that the input graph is already numbered from 0 to n-1 (improves performance)" << endl;
    cerr << "[-p|--partition] file\tstart the computation with a given partition instead of the trivial partition." << endl;
    cerr << "\tfile must contain lines \"node community\"." << endl;
    cerr << "[-q|--precision] eps\ta given pass stops when the modularity is increased by less than epsilon." << endl;
    cerr << "[-l|--level] k\tdisplays the graph of level k rather than the hierachical structure." << endl;
    cerr << "\tif k=-1 then displays the hierarchical structure rather than the graph at a given level (default)." << endl;
    cerr << "\tif k=-2 then displays the last level computed" << endl;
    cerr << "[-v|--verbose]\tverbose mode: gives computation time, information about the hierarchy and modularity." << endl;
    cerr << "[-h|--help]\tshow this usage message." << endl;
    exit(0);
}

void parse_args(int argc, char **argv) {
    if(argc == 1) {
        usage(argv[0]);
        exit(1);
    }

    static struct option long_options[] =
    {
        /* These options don’t set a flag.
           We distinguish them by their indices. */
        {"file",                required_argument, 0, 'f'},
        {"weighted",            no_argument,       0, 'w'},
        {"renumbered",          no_argument,       0, 'n'},
        {"reproducibility",     no_argument,       0, 'r'},
        {"partition",           required_argument, 0, 'p'},
        {"precision",           required_argument, 0, 'q'},
        {"level",               required_argument, 0, 'l'},
        {"verbose",             no_argument,       0, 'v'},
        {"help",                no_argument,       0, 'h'},
        {0, 0, 0, 0}
    }; 

    opterr = 0;
    // Parsing arguments using getopt 
    int arg;
    int option_index = 0;
    bool isnumber=true;
    // The first colon allows to separate error messages for missing argument and unknown options 
    while ((arg = getopt_long(argc, argv, ":f:rwnp:q:l:vh", long_options, &option_index)) != -1) {
        switch (arg) {
            case 0:
                break;
            case 'f':
                /* TODO: hack to handle flag recognized as argument, needs a better solution 
                 * for instance ./bin/community -f -a will assume that filename is "-a"
                 */
                if(optarg[0]=='-') {
                    cerr << "Option -f|--file requires an argument (input graph)" << endl;
                    exit(1);
                }
                else {
                    filename = optarg;
                    break;
                }
            case 'r':
                reproducibility = true;
                break;
            case 'w':
                weighted = true;
                break;
            case 'n':
                renumbering = false;
                break;
            case 'p':
                if(optarg[0]=='-') {
                    cerr << "Option -p|--partition requires an argument (partition file name)" << endl;
                    exit(1);
                }
                else {
                    filename_part = optarg;
                    break;
                }
            case 'q':
                if(optarg[0]=='-') {
                    cerr << "Option -q|--precision requires an argument (modularity gain threshold)" << endl;
                    exit(1);
                }
                else {
                    precision = atof(optarg);
                    break;
                }
            case 'l':
                isnumber = (string(optarg).find_first_not_of("0123456789") == string::npos);
                if(!isnumber && string(optarg)!="-1" && string(optarg)!="-2") {
                    cerr << "Option -l|--level expects a number as argument" << endl; 
                    exit(1);
                }
                else {
                    display_level = atoi(optarg);
                    break;
                }
            case 'v':
                verbose = true;
                break;
            case 'h':
                usage(argv[0]);
                break;
            case ':':
                if (optopt == 'f')
                    cerr << "Option -f|--file requires an argument (input graph)" << endl;
                else if (optopt == 'p')
                    cerr << "Option -p|--partition requires an argument (partition file name)" << endl;
                else if (optopt == 'q')
                    cerr << "Option -q|--precision requires an argument (modularity gain threshold)" << endl;
                else if (optopt == 'l')
                    cerr << "Option -l|--level requires an argument (level to display)" << endl;
                break;
            case '?':
                if (isprint (optopt))
                    cerr << "Unknown option '-" << (char)optopt << "'." << endl;
                else
                    cerr << "Unknown option character '\\x " << (char)optopt << "'." << endl;
                exit(1);
            default:
                abort();
        }
    }
}

#endif // UTILS_HPP
