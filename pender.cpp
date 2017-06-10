#include "pender.h"


bool Pender::alreadyProcessed(QString filename){
    QFile inputFile(filename);

    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);

        bool found_al_column = false;

        QString line ="##";

        while (!in.atEnd()){
            line = in.readLine();

            if (line.length() < 2) continue;

            // Skip any missed headers (just in case..)
            if (line[0]=='#') continue;

            // Test first data line
            QStringList tokens = line.split('\t');

            QStringList format_col = tokens.at(FORMAT_INDEX).split(':'); //"AL:GY:ZYG:ETC"

            if (format_col.contains(AL_ID)){
                found_al_column = true;
            }
            break;
        }
        inputFile.close();

        if (found_al_column){
            return true;
        }
    }
    return false;
}



void Pender::appendToVCF(QString colfile)
{
    QFile inputFile(colfile);
    uint numlines = FileUtils::countlines(colfile);
    uint countline = 0;

    QString id = FileUtils::extractID(colfile);

    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);
        QString line ="##";

        FileUtils::handleHeaders(line, in, countline, rej, out);

        //Begin parsing
        while (!in.atEnd()){
            countline ++;
            if(countline%2000==0){
                cerr << "\r\t" << id.toUtf8().data() << ":\t" << ((100*countline)/numlines) << "%   " << flush;
            }

            line = in.readLine();
            QStringList tokens = line.split('\t');

            //Begin processing
            QString chrom = tokens[0].trimmed();
            uint pos = QString(tokens[1]).trimmed().toUInt();


            bool foundThisGene=false;
            bool foundGene=false;
            bool foundIso=false;
            bool foundExtras=false;
            bool stop=false;

            QStringList allgeneKeys = chromemap[chrom].keys();
            GeneNameMap &allgenes = chromemap[chrom];

            QMap<QString,bool> genename; // out line

            for (int i=0; (i< allgeneKeys.length() && !stop); i++){
                QString &gene_name = allgeneKeys[i];
                GeneHolder *gh = allgenes[gene_name];

                // Variant is within gene bounds?
                if ((gh->maxmin.min <= pos) && (pos <= gh->maxmin.max)){
                    foundGene = true;
                    foundThisGene = true;
                    bool foundThisIso=false;

                    QStringList isogeneKeys = gh->isos.keys();
                    IsoformMap &isogenes = gh->isos;

                    // Process each isoform as it's own gene
                    for(int j=0; (j < isogeneKeys.length()) && !stop; j++){
                        QString gene = isogeneKeys[j];
                        IsoHolder *ish = isogenes[gene];

                        // variant within isoform bounds?
                        if ((ish->maxmin.min <= pos) && (pos <= ish->maxmin.max)){
                            foundIso = true;
                            foundThisIso = true;
                            bool foundThisExtra = false;

                            QStringList extraKeys = ish->extras.keys();
                            RegionMap &extras_ref = ish->extras;


                            for(int k=0; k < extraKeys.length(); k++){
                                QString extra = extraKeys[k];
                                MaxMin &mm = extras_ref[extra];

                                if ((mm.min <= pos) && (pos <= mm.max)){
                                   foundExtras = true;
                                   foundThisExtra = true;

                                   QString fullgene = gene + "|" + extra;
                                   // duplicates get overwritten
                                   genename[fullgene] = true;
                                }
                            }

                            // Only if input map has introns will they be included
                            /*if(!foundThisExtra){
                                QString fullgene = gene + "|Intron";
                                genename[fullgene] = true;
                            }*/
                        }
                    }
//                    if(!foundThisIso){
//                        gene.append("_NOTISOFORM_INTERGENIC");
//                        if(genename=="") genename = gene;
//                        else genename.append(gene);
//                    }
                }
            }

            //empty genelist = intergenic
            if (genename.size() == 0 ) genename[INTERGENIC] = true;

            // Build genename line
            QStringList genename_keys = genename.keys();
            QString out_genename = genename_keys.join(",");

            tokens[FORMAT_INDEX].append(':').append(AL_ID);
            tokens[INDIV_START_INDEX].append(':').append(out_genename);

            //Build line
            QString lineout;
            for (int k=0; k < tokens.length(); k++){
                lineout.append(tokens[k].trimmed()).append('\t');
            }
            lineout = lineout.trimmed();

            if(foundGene){
                if(foundIso){
                    if(foundExtras) out << lineout.toUtf8().data() << endl;
                    else ((this->keepall)?out:rej) << lineout.toUtf8().data() << endl;
                }
                else ((this->keepall)?out:rej) << lineout.toUtf8().data() << endl;
            } else ((this->keepall)?out:rej) << lineout.toUtf8().data() << endl;
        }
    }
    cerr << "\r\t" << id.toUtf8().data() << ":\t" << ((100*countline)/numlines) << "%   " << endl;

    inputFile.close();
    out.close();
    rej.close();
}



