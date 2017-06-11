#ifndef TYPES_H
#define TYPES_H

#include <QMap>

typedef unsigned int uint;

typedef struct {
    uint pos1;
    uint pos2;
} upair;

//    typedef struct {
//        upair coding;
//        upair transc;
//        bool direction;
//    } genedee;

class MaxMin{
public:
    uint max, min;
    MaxMin(){MaxMin(0,0);}
    MaxMin(uint mi, uint mx){min=mi;max=mx;}
    void updateMaxMin(uint x1, uint x2){
        min = (x1<min)?x1:min;
        max = (x2>max)?x2:max;
    }
};
typedef QMap<QString, MaxMin> RegionMap;


class IsoHolder{
public:
    MaxMin maxmin;
    RegionMap extras;
    IsoHolder(RegionMap extras, MaxMin maxmin){
        this->extras = extras;
        this->maxmin = maxmin;
    }
};
typedef QMap<QString, IsoHolder*> IsoformMap;

class GeneHolder{
public:
    MaxMin maxmin;
    IsoformMap isos;
    GeneHolder(IsoformMap genesAndIsos, MaxMin maxmin){
        this->isos = genesAndIsos;
        this->maxmin = maxmin;
    }
};
typedef QMap<QString, GeneHolder*> GeneNameMap;
typedef QMap<QString, GeneNameMap> ChromosomeMap ;

#endif
