#include "pender.h"


void Pender::appendToVCF()
{
    QFile inputFile(colfile);
    uint numlines = countlines(colfile);
    uint countline = 0;

    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);
        QString line ="##";

        handleHeaders(line, in, countline, rej);

        bool alreadychecked=false;

        //Begin parsing
        while (!in.atEnd()){
            countline ++;
            if(countline%2000==0)
                cerr << "\r\t\tColFile: " << ((100*countline)/numlines) << "%   " << flush;

            line = in.readLine();
            QStringList tokens = line.split('\t');

            //Check if file has already been processed
            if (!alreadychecked){
                QStringList format_col = tokens.at(FORMAT_INDEX).split(':'); //"AL:GY:ZYG:ETC"
                if (format_col.contains(AL_ID)){
                    cerr << colfile.toUtf8().data() << " has already been processed! Printing out everything." << endl;
                    cout << line.toUtf8().data() << endl;
                    while (!in.atEnd()){
                        cout << in.readLine().toUtf8().data() << endl;
                    }
                    exit(0);
                }
                alreadychecked = true;
            }
            //Begin processing

            QString chrom = tokens[0].trimmed();
            uint pos = QString(tokens[1]).trimmed().toUInt();

            QString genename="";

            bool foundThisGene=false;
            bool foundGene=false;
            bool foundIso=false;
            bool foundExtras=false;
            bool stop=false;

            QStringList allgeneKeys = chromemap[chrom].keys();
            QMap<QString, GeneHolder*> &allgenes = chromemap[chrom];

            for (int i=0; (i< allgeneKeys.length() && !stop); i++){
                QString &gene_name = allgeneKeys[i];
                GeneHolder *gh = allgenes[gene_name];

                if ((gh->maxmin.min <= pos) && (pos <= gh->maxmin.max)){
                    foundGene = true;
                    foundThisGene = true;
                    bool foundThisIso=false;

                    QStringList isogeneKeys = gh->isos.keys();
                    QMap<QString, IsoHolder*> &isogenes = gh->isos;

                    for(int j=0; (j < isogeneKeys.length()) && !stop; j++){
                        QString gene = isogeneKeys[j];
                        IsoHolder *ish = isogenes[gene];

                        if ((ish->maxmin.min <= pos) && (pos <= ish->maxmin.max)){
                            foundIso = true;
                            foundThisIso = true;
                            bool foundThisExtra = false;

                            QStringList extraKeys = ish->extras.keys();
                            QMap<QString, MaxMin> &extras = ish->extras;


                            for(int k=0; k < extraKeys.length(); k++){
                                QString extra = extraKeys[k];
                                MaxMin &mm = extras[extra];

                                if ((mm.min <= pos) && (pos <= mm.max)){
                                   foundExtras = true;
                                   foundThisExtra = true;

                                   gene.append("|"+extra);

                                   if(!appendmultiple){
                                       genename = gene;
                                       stop = true;
                                       break;
                                   }
                                   else {
                                       if(genename=="") genename = gene;
                                       else genename.append(","+gene);
                                   }
                                }
                            }
                            if(!foundThisExtra){
                                gene.append("|Intron");
                                if(genename=="") genename = gene;
                                else genename.append(","+gene);
                            }
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
            if (genename.length()==0) genename.append(INTERGENIC);

            tokens[FORMAT_INDEX].append(':').append(AL_ID);
            tokens[INDIV_START_INDEX].append(':').append(genename);

            //Build line
            QString lineout;
            for (int k=0; k < tokens.length(); k++){
                lineout.append(tokens[k].trimmed()).append('\t');
            }
            lineout = lineout.trimmed();

            if(foundGene){
                if(foundIso){
                    if(foundExtras) cout << lineout.toUtf8().data() << endl;
                    else ((this->keepall)?cout:rej) << lineout.toUtf8().data() << endl;
                }
                else ((this->keepall)?cout:rej) << lineout.toUtf8().data() << endl;
            } else ((this->keepall)?cout:rej) << lineout.toUtf8().data() << endl;
        }
    }
    cerr << "\r\t\tColFile: 100%         " << endl;
    inputFile.close();
}


void Pender::populateMap()
{
    QFile inputFile(mapfile);

    uint numlines = countlines(mapfile);
    uint count = 0;

    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);
        while(!in.atEnd()){
            ++count;
            //cerr << "\rMapFile: " << count << flush;
            if (count%300==0)
                cerr << "\rMapFile: " << ((100*count)/numlines) << '%' << flush;

            QStringList tokens = in.readLine().split('\t');

            QString chrom = tokens[0].trimmed();

            uint pos1 = QString(tokens[1]).trimmed().toUInt();
            uint pos2 = QString(tokens[2]).trimmed().toUInt();

            /*
            Possible configs:
            =================
            RPL23AP82|Exon1_5'UTR
            XKR3|Exon1
            PI4KA-ISOF1|Exon13
            RABL2B-ISOF4|3'UTR
            SHANK3
            LOC101101776|Intron1
            OR11H1-CCT8L2|Intergenic
            BCL2L13-ISOF10|Exon2_SpliceD
            */
            QString genxon = QString(tokens[3]).trimmed();
            QString gene = genxon;  // Assuming no '|'
            QString geneNoIso = gene.split("-ISO")[0].split('|')[0].trimmed();

            QString extra= ""; // Nothing extra             Exon_Splice,Intron,Intergenic

            if(genxon.contains('|')){
                QStringList g = genxon.split('|');

                gene = g[0].trimmed(); //Update to real gene
                extra = (g.length() > 1)?(g[1].trimmed()):"";
            }

            if (skipbadscore && ((tokens.length()>4))){
                QString score = QString(tokens[4]).trimmed();
                if (score!="cC") continue;
            }

            if (!chromemap.contains(chrom)){
                QMap<QString,GeneHolder* > genemap;
                chromemap[chrom] = genemap;
            }
            else {
                //Chrom exists, check gene names
                QMap<QString, GeneHolder*> &allgenenames = chromemap[chrom];

                // All genes in that chromosome -- names
                if (!allgenenames.contains(geneNoIso)){
                    QMap<QString, IsoHolder*> isos;
                    GeneHolder *gh = new GeneHolder(isos, MaxMin(pos1, pos2));
                    allgenenames[geneNoIso] = gh;

                }
                else { // All genes within that gene -- isoforms and all
                    GeneHolder *isos_gene = allgenenames[geneNoIso];
                    isos_gene->maxmin.updateMaxMin(pos1,pos2);
                    QMap<QString, IsoHolder*> &isos = isos_gene->isos;

                    if(!isos.contains(gene)){
                        QMap<QString,MaxMin> extras;
                        isos[gene] = new IsoHolder(extras, MaxMin(pos1,pos2));
                    }
                    else { //All extras within that isoform
                        IsoHolder *isogenes = isos[gene];
                        isogenes->maxmin.updateMaxMin(pos1,pos2);
                        QMap<QString,MaxMin> &extras = isogenes->extras;

                        if(extra==""){
                            if(!extras.contains(extra)){
                                extras[extra] = MaxMin(isogenes->maxmin.min, isogenes->maxmin.max);
                            }
                        } else {
                            if(!extras.contains(extra)){
                                extras[extra] = MaxMin(pos1,pos2);
                            }
                        }
                    }
                }
            }
        }
    }
    cerr << "\rMapFile: 100% \t" << flush;
//    QStringList keys1 = chromemap.keys();
//    cerr << "Number of chroms:" << keys1.length() << endl;
//    QStringList keys2 = chromemap[keys1[0]].keys();
//    cerr << "Number of genes in " << keys1[0].toUtf8().data() << ": " << keys2.length() << endl;
//    QStringList keys3 = chromemap[keys1[0]][keys2[3]]->isos.keys();
//    cerr << "Number of isoforms in " << keys2[3].toUtf8().data() << ": " << keys3.length() << endl;
//    cerr << keys2[3].toUtf8().data() << "  "
//            << chromemap[keys1[0]][keys2[3]]->maxmin.min << "__"
//            << chromemap[keys1[0]][keys2[3]]->maxmin.max << endl;
//    QStringList keys4 = chromemap[keys1[0]][keys2[3]]->isos[keys3[0]]->extras.keys();
//    cerr << "Number of extras in " << keys3[0].toUtf8().data() << ": " << keys4.length() << endl;
    inputFile.close();
}


