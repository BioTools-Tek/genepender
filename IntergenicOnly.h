#include "MapInterface.h"

#ifndef INTERGENICONLY_H
#define INTERGENICONLY_H

//#define DEBUG_MAP


class IntergenicOnly : public MapInterface {
public:
    IntergenicOnly(QString &map, QString &colf, QString &rejects, bool &append, bool &skipbad){
        mapfile = map;
        colfile = colf;
        appendmultiple = append;
        skipbadscore = skipbad;
        //Open Channel
        rej.open(rejects.toUtf8().data());

        populateMap();
        appendToVCF();
    }

    ~IntergenicOnly(){ rej.close();}

    ChromosomeMap chromemap;

    void appendToVCF();
    void populateMap();
};


void IntergenicOnly::appendToVCF(){
    QFile inputFile(colfile);

    uint numlines = countlines(colfile);
    uint countline = 0;

    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);
        QString line ="##";

        handleHeaders(line, in, countline, rej, out);

        //Begin parsing
        while (!in.atEnd()){
            countline ++;
            if(countline%100==0)
                cerr << "\r\t\tColFile: " << ((100*countline)/numlines) << "%   " << flush;

            line = in.readLine();
            QStringList tokens = line.split('\t');

            QString chrom = tokens[0].trimmed();
            uint pos = QString(tokens[1]).trimmed().toUInt();

            bool foundgene = false;
            QString genename="";

            QStringList genes = chromemap[chrom].keys();

            for (int i=0; i< genes.length(); i++){
                QString &gene = genes[i];
                genedee &datar = chromemap[chrom][gene];
                upair &positions = datar.coding;
                upair &transcrip = datar.transc;

                if (( transcrip.pos1 <= pos ) && (pos <= transcrip.pos2 )){
                    //Within transcription region.... good.
                    bool utr = false;

                    //StartUTR
                    if ( pos < positions.pos1 ){
                        utr = true;
                        if(datar.direction) gene.append(UTR5);
                        else gene.append(UTR3);
                    }
                    else if (pos > positions.pos2) {
                        utr = true;
                        if(datar.direction) gene.append(UTR3);
                        else gene.append(UTR5);
                    }

//                    if(!utr) foundgene = true;  //Print UTR's to rejects file (ONLY IF ALL APPENDED GENES ARE UTR. IF JUST ONE ISNT --> COUT)
                    foundgene = true; //Keep UTR in output

                    if (!appendmultiple) { genename=gene; break; }
                    //Append
                    else {
                        if (genename=="") genename = gene;
                        else genename.append(","+gene);
                    }
                }
            }
            //No gene found for position, integenic
            if(!foundgene) genename.append(INTERGENIC);

            tokens[FORMAT_INDEX].append(':').append(GEL_ID);
            tokens[INDIV_START_INDEX].append(':').append(genename);


            //Build line
            QString lineout;
            for (int k=0; k < tokens.length(); k++){
                lineout.append(tokens[k].trimmed()).append('\t');
            }

            if (foundgene) out << lineout.toUtf8().data() << endl;
            else rej << lineout.toUtf8().data() << endl;
        }
    }
    cerr << "\r\t\tColFile: 100%         " << endl;
    inputFile.close();
}


void IntergenicOnly::populateMap()
{
    QFile inputFile(mapfile);

    uint numlines = countlines(mapfile);
    uint count = 0;

    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);
        while(!in.atEnd()){
            ++count;
            if (count%300==0)
                cerr << "\rMapFile: " << ((100*count)/numlines) << '%' << flush;

            QStringList tokens = in.readLine().split('\t');

            QString chrom = tokens[0];
            uint pos1 = QString(tokens[1]).trimmed().toUInt();
            uint pos2 = QString(tokens[2]).trimmed().toUInt();
            QString genxon = QString(tokens[3]);
            QString gene = genxon.split('_')[0].trimmed();

            if (skipbadscore && ((tokens.length()>4))){
                QString score = QString(tokens[4]).trimmed();
                if (score!="") continue;
            }

            if (!chromemap.contains(chrom)){
                QMap<QString,genedee > genemap;
                chromemap[chrom] = genemap;
            }
            else {
                //Chrom exists, check genemap
                QMap<QString, genedee > &genemap = chromemap[chrom];

                if (!genemap.contains(gene)){
                    //Gene doesn't exist in map
                    genedee datar;

                    upair code_positions;
                    code_positions.pos1 = pos1;
                    code_positions.pos2 = pos2;

                    //find txStart+txStop just once.
                    QStringList txes = tokens[6].trimmed().split('_');
                    uint txS = txes[0].toUInt();
                    uint txE = txes[1].toUInt();

                    upair trans_pos;
                    trans_pos.pos1 = txS;
                    trans_pos.pos2 = txE;

                    //Assign to datar
                    datar.direction = (tokens[5]=="+");
                    datar.coding = code_positions;
                    datar.transc = trans_pos;

                    genemap[gene] = datar;
                }
                else {
                    //Gene exists, update max min
                    genedee &datar = genemap[gene];
                    upair &code_positions = datar.coding;
                    if (pos1 < code_positions.pos1) code_positions.pos1 = pos1;
                    if (pos2 > code_positions.pos2) code_positions.pos2 = pos2;
                }
            }
        }
    }
    cerr << "\rMapFile: 100% " << flush;

#ifdef DEBUG_MAP
    QStringList chroms = chromemap.keys();
    for (int ll=0; ll< chroms.length(); ll ++){
        cerr << "MAPSIZE " << chroms.at(ll).toUtf8().data() << ':' << chromemap[chroms[ll]].size() << "genes" << endl;
    }
#endif


    inputFile.close();
}


#endif // INTERGENICONLY_H
