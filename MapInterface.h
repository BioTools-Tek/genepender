#ifndef MAPINTERFACE_H
#define MAPINTERFACE_H

//STRING CONSTANTS
#define SPLICE "_SPLICE"
#define EXON "_EXON"
#define INTRON "_INTRON"
#define INTERGENIC "INTERGENIC"
#define UTR5 "_5UTR"
#define UTR3 "_3UTR"

#define AL_ID "AL"
#define HFORMAT "##FORMAT="
#define FILLER ",Number=.,Type=String,Description="
#define HEADER_ALIST HFORMAT "<ID=" AL_ID
#define HEADER_ALIST_FULL HEADER_ALIST FILLER "\"List of Genes with Exons\">"


#include <QString>
#include <QMap>
#include <iostream>
#include <fstream>
#include <QTextStream>
#include <QStringList>
#include <QFile>

using namespace std;

static int FORMAT_INDEX = 8; // Default
static int INDIV_START_INDEX = 9; // Default

class MapInterface {
public:

    typedef struct {
        uint pos1;
        uint pos2;
    } upair;

//    typedef struct {
//        upair coding;
//        upair transc;
//        bool direction;
//    } genedee;


    QString mapfile;
    QString colfile;
    QString outfile;
    ofstream rej;
    ofstream out;
    bool appendmultiple;
    bool skipbadscore;
    bool keepall;
    bool forceprocs;

    static uint countlines(QString file){
        uint lines=0;
        //For lines in infile
        QFile inputFile(file);
        if (inputFile.open(QIODevice::ReadOnly))
        {
            while (!inputFile.atEnd()){
                inputFile.readLine();
                lines ++;
            }
        }
        return lines;
    }

    static bool handleHeaders(QString &line, QTextStream &in, uint &countline,
                              ofstream &rej, ofstream &out)
    {
        line="##";

        qint64 bufferpos = 0;

        bool found_header_alist = false;
        short found_format = -1;

        //Print headers
        do{
            line = in.readLine();
            if (line[0]!='#') break;

            //Look for ##
            if (line.startsWith(HFORMAT)){
                found_format = 0; //Found
                if (line.startsWith(HEADER_ALIST)) found_header_alist = true;
            }
            else if (line.startsWith("#CHROM")){
                QStringList tokes = line.split('\t');
                bool form = false;
                for (int p=0; p < tokes.length(); p++){
                    if (tokes[p].toUpper().contains("FORMAT")) {
                        FORMAT_INDEX = p;
                        INDIV_START_INDEX = p+1;
                        form = true;
                        break;
                    }
                }
                if (!form) cerr << "Could not find FORMAT column! Assuming column:" << (FORMAT_INDEX + 1) << endl;

            }


            else {
                if(found_format==0){ // Found format but now it's ended
                    if (!found_header_alist) out << HEADER_ALIST_FULL << endl; //Never found header, print new one
                    found_format=-2; // so that we know that it at least exists
                }
            }
            out << line.toUtf8().data() << endl; // Output to rejects and out too
            rej <<  line.toUtf8().data() << endl;

            countline++;
            bufferpos = in.pos();
        } while(!in.atEnd());

        //Never found ##FORMAT in header
        if(found_format==-1){
            cerr << "No Format line in header!\nPrinting new one anyway..." << endl;
            if (!found_header_alist) out << HEADER_ALIST_FULL << endl; //Never found header, print new one
        }
        //Go to end of headers
        in.seek(bufferpos);

        return found_header_alist;
    }
};


#endif // MAPINTERFACE_H
