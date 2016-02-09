//#include "IntergenicOnly.h"
//#include "InterAndExonic.h"*/

// IntergenicOnly and InterAndExonic are DEPRECATED

#include "pender.h"

//Flags
#define MULTIPLE_ID "--multiple"
#define SCORES_ID "--allscores"
#define KEEPINT_ID "--keepall"
#define FORCE_ID "--force"

void usage(){
    cerr << "v2016_02_09" << endl;
    cerr << "Appends a new column containing gene names, whilst filtering out intergenic regions (i.e. where there are no genes), and can also optionally filter out introns.\n" << endl;
    cerr << "Usage: genepender <genemap> <vcf in> <rejects> <vcf out> <[OPTIONS]" << endl;
    cerr << endl;
    cerr << "where a <column file> is any file which has a <chr>[TAB]<index>[TAB]<etc...> layout" << endl;
    cerr << "and <rejects> is a file created to store intergenic and intron lines from column file" << endl;
    cerr << endl;
    cerr << "OPTIONS:" << endl;
    cerr << KEEPINT_ID << "\tAlso keeps variants that only fall within intron/intergenic regions." << endl;
    cerr << MULTIPLE_ID << "\tPrints all gene names if index bisects multiple genes." << endl;
    cerr << FORCE_ID << "\t\tForces reprocessing even if VCF header present in either in/out file" << endl;
    cerr << SCORES_ID << "\tDont skip exons with poor map scores" << endl;
    cerr <<            "\t\tTypically: uU, iI, nN, and wW" << endl;
    cerr << endl;
    exit(-1);
}

int main(int argc, char *argv[])
{
    bool    appendmultiple = false,
            skipbad = true,
            keepints = false,
            force = false;

    if (argc<5) usage();
    else if (argc > 6){
        for (int i=6; i < argc; i++){
            QString arg = argv[i];
            if ( arg==MULTIPLE_ID) appendmultiple = true;
            else if ( arg==SCORES_ID) skipbad = false;
            else if ( arg==FORCE_ID) force = true;
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
    QString outfile = argv[4];

    Pender(mapfile, colfile, rejects, outfile,
           appendmultiple, skipbad, keepints, force);

    //Eh, messy. Have to pass bool to static handleHeaders
    if (keepints){
        QString comm= "rm "+rejects;
        system(comm.toUtf8().data());
    }
    return 0;
}

