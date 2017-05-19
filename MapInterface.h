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

#endif // MAPINTERFACE_H
