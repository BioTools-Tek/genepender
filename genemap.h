#include "pender.h"

class GeneMap {
public:
    ChromosomeMap chromemap;
    //chrom-->Gene name --> maxmin, and geneSo

    GeneMap(QString filename){
        this->populateMap(filename);
    }

private:
    void populateMap(QString mapfile){
        QFile inputFile(mapfile);

        uint numlines = FileUtils::countlines(mapfile);
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

                // depreciated
                /*if (skipbadscore && ((tokens.length()>4))){
                    QString score = QString(tokens[4]).trimmed();
                    if (score!="cC") continue;
                }*/

                if (!chromemap.contains(chrom)){
                    GeneNameMap genemap;
                    chromemap[chrom] = genemap;
                }
                else {
                    //Chrom exists, check gene names
                    GeneNameMap &allgenenames = chromemap[chrom];

                    // All genes in that chromosome -- names
                    if (!allgenenames.contains(geneNoIso)){
                        IsoformMap isos;
                        GeneHolder *gh = new GeneHolder(isos, MaxMin(pos1, pos2));
                        allgenenames[geneNoIso] = gh;

                    }
                    else { // All genes within that gene -- isoforms and all
                        GeneHolder *isos_gene = allgenenames[geneNoIso];
                        isos_gene->maxmin.updateMaxMin(pos1,pos2);
                        IsoformMap &isos = isos_gene->isos;

                        if(!isos.contains(gene)){
                            RegionMap extras;
                            isos[gene] = new IsoHolder(extras, MaxMin(pos1,pos2));
                        }
                        else { //All extras within that isoform
                            IsoHolder *isogenes = isos[gene];
                            isogenes->maxmin.updateMaxMin(pos1,pos2);
                            RegionMap &extra_refs = isogenes->extras;

                            if(extra==""){
                                if(!extra_refs.contains(extra)){
                                    extra_refs[extra] = MaxMin(isogenes->maxmin.min, isogenes->maxmin.max);
                                }
                            } else {
                                if(!extra_refs.contains(extra)){
                                    extra_refs[extra] = MaxMin(pos1,pos2);
                                }
                            }
                        }
                    }
                }
            }
        }
        cerr << "\rMapFile: 100% \n" << flush;
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
};
