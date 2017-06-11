#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <QFile>
#include <QDir>
#include <QTextStream>
#include <QStringList>

static int FORMAT_INDEX = 8; // Default
static int INDIV_START_INDEX = 9; // Default

class FileUtils {
public:
    static bool makeFolder(QString dirname){
        //cerr << "Request path:" << dirname.toUtf8().data() << endl;
        QDir dir(dirname);
        if (!dir.exists()){
            //cerr << "Making path:" << dir.dirName().toUtf8().data() << endl;
            return dir.mkpath(".");
        }
        return false;
    }

    static QString makePrefixRootDir(QString input_filename, QString output_dir, QString prefix)
    {
        int ind = input_filename.lastIndexOf(prefix);
        if (ind == -1){
            cerr << "Prefix not in filename: " << prefix.toUtf8().data()
                 << " " << input_filename.toUtf8().data() << endl;
            exit(-1);
        }

        QString family_structure_file = input_filename.mid(ind);
        QStringList family_structure_dir = family_structure_file.split("/");
        family_structure_dir.removeLast();
        QString family_structure = family_structure_dir.join("/");

        QString new_output_dir = QDir(output_dir).filePath(family_structure);

        bool res = FileUtils::makeFolder(new_output_dir);
        //cerr << "NEW OUT = " << new_output_dir.toUtf8().data() << "-- " << res << endl;

        return new_output_dir;
    }

    static QString extractID(QString filename){
        return filename.split("/").last().split(".").first();
    }

    static QString suffixToFilename(QString filename, QString suff){
        QStringList f_split = filename.split('.') ;
        f_split.insert( f_split.length()-1, suff ); // 25.blah.space.vcf --> 25.blah.space.suff.vcf
        return f_split.join(".");
    }

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




#endif
