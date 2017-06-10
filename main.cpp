//#include "IntergenicOnly.h"
//#include "InterAndExonic.h"*/

// IntergenicOnly and InterAndExonic are DEPRECATED

#include "genemap.h"

//Flags
//#define MULTIPLE_ID "--multiple" -- default
//#define SCORES_ID "--allscores" -- default
#define KEEPINT_ID "--keepall"
#define FORCE_ID "--force"
#define OUTPATH_ID "--output-folder="
#define PREFEXTR_ID "--prefix-extract="

void usage(){
    cerr << "2.6v20170522" << endl;
    cerr << endl;
    cerr << "Appends a new column containing gene names, whilst filtering out intergenic regions (i.e. where there are no genes), and can also optionally filter out introns.\n" << endl;
    cerr << "Usage: genepender <genemap> <[VCFS]> <[OPTIONS]" << endl;
    cerr << endl;
    cerr << "where a <column file> is any file which has a <chr>[TAB]<index>[TAB]<etc...> layout." << endl;
    cerr << "Two files are created: one with a .genes.rejects.vcf suffix, and the other with a .genes.vcf suffix at the source directory" << endl;
    cerr << endl;
    cerr << "OPTIONS:" << endl;
    cerr << KEEPINT_ID << "\tAlso keeps variants that only fall within intron/intergenic regions." << endl;
    cerr << FORCE_ID << "\t\tForces reprocessing even if VCF header present in either in/out file" << endl;
    cerr << OUTPATH_ID << "<name>\tRoot output folder, source folder otherwise" << endl;
    cerr << PREFEXTR_ID << "<name>\tPreserves directory structure starting from prefix value. Dumped into root of output folder otherwise." << endl;
    cerr << endl;
    exit(-1);
}

int main(int argc, char *argv[])
{
    bool keepints = false, force = false;
    QString outdir = "";
    QString prefext = "";

    QStringList files;

    for (int i=2; i < argc; i++){
        QString arg = argv[i];
        if      ( arg==FORCE_ID    ) force = true;
        else if ( arg==KEEPINT_ID  ) keepints = true;
        else if ( arg.startsWith(OUTPATH_ID)  ){ outdir = arg.split(OUTPATH_ID).last();}
        else if ( arg.startsWith(PREFEXTR_ID) ){ prefext= arg.split(PREFEXTR_ID).last();}
        else if ( arg.startsWith("--")){
            cerr << "Unable to parse:" << arg.toUtf8().data() << endl;
            exit(-1);
        }
        else {
            files.append(arg);
        }
    }
    if (files.length()==0) usage();

    QString mapfile = argv[1];
    GeneMap *gmp = new GeneMap(mapfile);

    for (int f=0; f < files.length(); f++)
    {
        Pender(gmp->chromemap, files[f], outdir, prefext, keepints, force);
    }

    return 0;
}

