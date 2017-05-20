#include "MapInterface.h"



#ifndef PENDER_H
#define PENDER_H

#include "FileUtils.h"
#include "Types.h"

class Pender {
public:
    ofstream rej;
    ofstream out;
    bool keepall;
    bool forceprocs;

    //Pender(gmp->chromemap, colfile, keepints, force);
    Pender(ChromosomeMap &chromemap, QString &colf, QString &outdir, QString &prefext, bool &keepints, bool &force_processing)
    {
        QString outf = FileUtils::suffixToFilename(colf, "genes");
        QString rejects = FileUtils::suffixToFilename(colf, "genes.rejects");

        QString new_outdir = outdir;

        if (prefext != ""){
            new_outdir = FileUtils::makePrefixRootDir(colf, outdir, prefext);
        }

        if (outdir == ""){
            new_outdir = "./";
        }

        outf = new_outdir + "/" + outf.split("/").last();
        rejects = new_outdir + "/" + rejects.split("/").last();


        //cerr << "outf = " << outf.toUtf8().data() << endl;
        //cerr << "rejects = " << rejects.toUtf8().data() << endl;

        keepall = keepints;
        forceprocs = force_processing;

        this->chromemap = chromemap;

        if (!forceprocs){
            bool processedInput  = alreadyProcessed(colf);
            bool processedOutput = alreadyProcessed(outf);

            if (processedInput){
                cerr << "Gene headers present in input. Skipping" << endl;
                exit(1);
            }
            if (processedOutput){
                cerr << "Gene headers present in output. Skipping." << endl;
                exit(2);
            }
        }

        // Either exitted due to previous processing, or go ahead given.

        //Open Channels
        rej.open(rejects.toUtf8().data());
        out.open(outf.toUtf8().data());

        // cerr << "Rejects: " << rejects.toUtf8().data() << endl;
        // cerr << "Outfile: " << outf.toUtf8().data() << endl;

        appendToVCF(colf);
    }



    ~Pender(){}

    ChromosomeMap chromemap;

    void appendToVCF(QString filename);
    void populateMap();
    bool alreadyProcessed(QString filename);
};


#endif // PENDER_H
