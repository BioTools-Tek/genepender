#include "MapInterface.h"

#ifndef INTERANDEXONIC_H
#define INTERANDEXONIC_H



class InterAndExonic : public MapInterface {
private:
    typedef struct {
        int utr_num;
        bool direction;
        upair maxmin;
        upair txStarSto;
        QMap<uint, upair > exons;
    } GeneHolder;
public:
    int splice_sites;

    InterAndExonic(QString &map, QString &colf, QString &rejects, bool &append, bool &skipbad, int &splice){
        mapfile = map;
        colfile = colf;
        appendmultiple = append;
        skipbadscore = skipbad;
        splice_sites = splice;
        //Open Channel
        rej.open(rejects.toUtf8().data());

        populateMap();
        appendToVCF();
    }

    ~InterAndExonic(){ rej.close(); }

//    [chrom, [gene, [  [max,min], [exon, [ pos1, pos2 ] ] ] ]]        [Qstring, [QString, [GeneHolder]]]
    QMap<QString, QMap<QString, GeneHolder> > chromemap;

    void appendToVCF();
    void populateMap();
};

void InterAndExonic::appendToVCF(){
    QFile inputFile(colfile);

    uint numlines = countlines(colfile);
    uint countline = 0;

    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);
        QString line="##";

        handleHeaders(line, in, countline, rej, true); // true print exon headers, find FORMAT_INDEX

        //Begin parsing
        while (!in.atEnd()){
            countline ++;
            if(countline%100==0)
                cerr << "\r\t\tColFile: " << ((100*countline)/numlines) << "%   " << flush;

            line = in.readLine();
            QStringList tokens = line.split('\t');

            QString chrom = tokens[0].trimmed();
            uint pos = QString(tokens[1]).trimmed().toUInt();

            QString genelist = "";
            QStringList genes = chromemap[chrom].keys();

            //Find single Gene+Exon
            bool foundOneExon=false;
            bool foundOneGene=false;

            //Find if pos is in Max-Min positions of gene
            for (int i=0; i< genes.length(); i++)
            {
                QString &gene = genes[i];
                GeneHolder &gh = chromemap[chrom][gene];

//                uint minsplice = gh.maxmin.pos1 - splice_sites;
//                uint maxsplice = gh.maxmin.pos2 + splice_sites;

                //Within transcript region?
                if (( gh.txStarSto.pos1 <= pos) && ( pos <= gh.txStarSto.pos1 )){
                    bool utr = false;

                    if ( pos < gh.maxmin.pos1 ){
                        if(gh.direction) gene.append(UTR5);
                        else gene.append(UTR3);
                        utr = true;
                    }
                    if ( pos > gh.maxmin.pos2){
                        if(gh.direction) gene.append(UTR3);
                        else gene.append(UTR5);
                        utr = true;
                    }
                    if(utr){
                        if(!appendmultiple) {genelist= gene; break;}
                        else{
                            if(genelist=="") genelist = gene;
                            else genelist.append(","+gene);
                        }
                    }

                    if(!utr){  // within actual gene
                        foundOneGene=true;
                        QList<uint> exons = gh.exons.keys();

                        bool foundThisExon = false;

                        for (int j=0; j < exons.length(); j++){
                            uint exon_number = exons[j];
                            upair &positions = gh.exons[exon_number];

                            //Found exon
                            uint minrange = positions.pos1 - splice_sites;
                            uint maxrange = positions.pos2 + splice_sites;

                            if ( (minrange <= pos) && (pos <= maxrange) ){
                                foundOneExon=true;
                                foundThisExon=true;

                                gene.append(EXON).append(QString::number(exon_number));

                                //Is it on the outskirts of the actual region? splice.
                                if (( pos < positions.pos1 ) || (pos > positions.pos2) ) gene.append(SPLICE).append(QString::number(splice_sites));

                                if (genelist=="") genelist = gene;
                                else genelist.append(","+gene);

    //                          A single position wont give multiple exons of the same gene - break.
                                break;
                            }
                        }

                        if (!foundThisExon) {
                            gene.append(INTRON);// This is being printed t oo
                            if (genelist=="") genelist = gene;
                            else genelist.append(","+gene);
                        }
                        else if (foundThisExon){
                            if(!appendmultiple) break;  // Found single exon, for single gene, job done.
                        }
                    }
                }
            }
            //Swap in new data
//            QString format = tokens[FORMAT_INDEX].trimmed(); //Append format
//            format.append(':').append(EXL_ID);

            //Empty genelist = Intergenic
            if (genelist.length()==0) genelist.append(INTERGENIC);

            tokens[FORMAT_INDEX].append(':').append(EXL_ID);
            tokens[INDIV_START_INDEX].append(':').append(genelist);

            //Build line
            QString lineout;
            for (int k=0; k < tokens.length(); k++){
                lineout.append(tokens[k].trimmed()).append('\t');
            }

            if (foundOneGene){
                if (foundOneExon){
                    //GeneName + Exon
                    cout << lineout.toUtf8().data() << endl;
                }
                else { //GeneName + Intron
                    rej << lineout.toUtf8().data() << endl;
                }
            }
            else { //NO genename
                //rej << line.toUtf8().data() << INTERGENIC << endl;
                rej << lineout.toUtf8().data() << endl;
            }
        }
    }
    cerr << "\r\t\tColFile: 100%         " << endl;
    inputFile.close();
}


void InterAndExonic::populateMap(){
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

            QString chrom = tokens[0];
            uint pos1 = QString(tokens[1]).trimmed().toUInt();
            uint pos2 = QString(tokens[2]).trimmed().toUInt();

            QString genxon = QString(tokens[3]);
            QStringList g = genxon.trimmed().split('_');           
            uint exon = g.last().split("Exon")[1].toUInt();

            int utr=-1;
            if (g.size()==3) utr = g[1].split("UTR")[0].toInt();

            QString gene = g[0].trimmed();

            if (skipbadscore && ((tokens.length()>4))){
                QString score = QString(tokens[4]).trimmed();
                if (score!="") continue;
            }

            if (!chromemap.contains(chrom)){
                QMap<QString,GeneHolder > genemap;
                chromemap[chrom] = genemap;
            }

            else {
                //Chrom exists, check genemap
                QMap<QString, GeneHolder > &genemap = chromemap[chrom];

                if (!genemap.contains(gene)){
                    GeneHolder gg;
                    gg.maxmin.pos1 = pos1;
                    gg.maxmin.pos2 = pos2;

                    //find txStart+txStop just once.
                    QStringList txes = tokens[6].trimmed().split('_');
                    gg.txStarSto.pos1 = txes[0].toUInt();
                    gg.txStarSto.pos2 = txes[1].toUInt();
                    genemap[gene] = gg;
                }

                //Gene exists, check exon
                else {
                    GeneHolder &hh = genemap[gene];
                    if (pos1 < hh.maxmin.pos1) hh.maxmin.pos1 = pos1;
                    if (pos2 > hh.maxmin.pos2) hh.maxmin.pos2 = pos2;

                    //Check if exon exists? Unlikely
                    //If so, quicker to overwrite duplicate data than to check for it
                    upair posse; posse.pos1 = pos1; posse.pos2 = pos2;
                    hh.exons[exon] = posse;
                }
            }
        }
    }
    cerr << "\rMapFile: 100% " << flush;
    inputFile.close();
}


#endif // INTERANDEXONIC_H
