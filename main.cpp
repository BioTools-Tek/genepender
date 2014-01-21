//#include "IntergenicOnly.h"
//#include "InterAndExonic.h"*/

// IntergenicOnly and InterAndExonic are DEPRECATED

#include "pender.h"

//Flags
#define MULTIPLE_ID "--multiple"
#define SCORES_ID "--allscores"
#define KEEPINT_ID "--keepall"

void usage(){
    cerr << "v2013_07_16" << endl;
    cerr << "Appends a new column containing gene names, whilst filtering out intergenic regions (i.e. where there are no genes), and can also optionally filter out introns.\n" << endl;
    cerr << "Usage: genepender <genemap> <column file> <rejects> [OPTIONS]" << endl;
    cerr << endl;
    cerr << "where a <column file> is any file which has a <chr>[TAB]<index>[TAB]<etc...> layout" << endl;
    cerr << "and <rejects> is a file created to store intergenic and intron lines from column file" << endl;
    cerr << endl;
    cerr << "OPTIONS:" << endl;
    cerr << KEEPINT_ID << "\tKeeps variants that fall ONLY in introns and/or intergenic regions too" << endl;
    cerr << MULTIPLE_ID << "\tPrints all gene names if index bisects multiple genes." << endl;
    cerr << SCORES_ID << "\tDont skip exons with poor map scores" << endl;
    cerr <<            "\t\tTypically: uU, iI, nN, and wW" << endl;
    cerr << endl;
    exit(-1);
}

int main(int argc, char *argv[])
{
    bool appendmultiple = false;
    bool skipbad = true;
    bool keepints = false;

    if (argc<4) usage();
    else if (argc >4){
        for (int i=4; i < argc; i++){
            QString arg = argv[i];
            if ( arg==MULTIPLE_ID) appendmultiple = true;
            else if ( arg==SCORES_ID) skipbad = false;
            else if ( arg==KEEPINT_ID) keepints = true;
            else {
                cerr << "Unable to parse:" << arg.toUtf8().data() << endl;
                exit(-1);
            }
        }
    }

    QString mapfile = argv[1];
    QString colfile = argv[2];
    QString rejects = argv[3];

    Pender(mapfile, colfile, rejects,
           appendmultiple, skipbad, keepints);

    //Eh, messy. Have to pass bool to static handleHeaders
    if (keepints){
        QString comm= "rm "+rejects;
        system(comm.toUtf8().data());
    }
    return 0;
}

