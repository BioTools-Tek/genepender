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

                QStringList tokens = in.readLine().split('\t');

                QString chrom = tokens[0].trimmed();

                uint pos1 = QString(tokens[1]).trimmed().toUInt();
                uint pos2 = QString(tokens[2]).trimmed().toUInt();

                QString genxon = QString(tokens[3]).trimmed();
                QString gene = genxon;  // Assuming no '|'
                QString geneNoIso = gene.split("-ISO")[0].split('|')[0].trimmed();

                QString extra= ""; // Nothing extra             Exon_Splice,Intron,Intergenic

                /*
                Possible configs:
                =================
                Gene1
                Gene1|Exon1
                Gene1-ISOF1|Exon13
                Gene1|Exon1_5'UTR
                Gene1|Intron1
                Gene1-ISOF4|3'UTR
                Gene1-ISOF10|Exon2_SpliceD
                Gene1-Gene2|Intergenic
                */

                if (genxon.contains('|')){
                    QStringList g = genxon.split('|');

                    gene = g[0].trimmed(); //Update to real gene
                    extra = (g.length() > 1)?(g[1].trimmed()):"";
                }

                // New chrom, new map
                if (!chromemap.contains(chrom)){
                    GeneNameMap genemap;
                    chromemap[chrom] = genemap;
                }

                //Chrom now exists, check gene names
                GeneNameMap &allgenenames = chromemap[chrom];

                // New gene, new geneholder in chrom, insert current as Max Min bounds of gene
                if (!allgenenames.contains(geneNoIso)){
                    IsoformMap isos;
                    GeneHolder *gh = new GeneHolder(isos, MaxMin(pos1, pos2));
                    allgenenames[geneNoIso] = gh;
                }


                // All genes within that gene -- isoforms and all
                GeneHolder *isos_gene = allgenenames[geneNoIso];
                isos_gene->maxmin.updateMaxMin(pos1,pos2);
                IsoformMap &isos = isos_gene->isos;

                // Is gene(isoform) in isoform map?
                if(!isos.contains(gene)){
                    RegionMap extras;
                    isos[gene] = new IsoHolder(extras, MaxMin(pos1,pos2));
                }


                //All extras within that isoform
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
                if (++count%300==0) cerr << "\rMapFile: " << ((100*count)/numlines) << '%' << flush;
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
